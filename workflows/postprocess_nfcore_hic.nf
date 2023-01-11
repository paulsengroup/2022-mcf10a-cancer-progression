#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    input_dirs = Channel.fromPath(params.nfcore_stage_dirs, type: 'dir')
    
    sample_names = [:]
    conditions = [:]

    input_dirs.subscribe onNext: {

        sample_names[it.toString()] = file(it).getName()
        conditions[it.toString()] = file(it).getName().replaceAll(/HiC_\d+_(.*)$/, '$1')
    }
    input_dirs_ = input_dirs.collect()

    input_coolers = input_dirs.flatten().map {
                                        tuple(conditions[it.toString()],
                                              sample_names[it.toString()],
                                              file("${it}/contact_maps/raw/cool/*.cool", type: "file", checkIfExists: true))
                                     }
    input_bwt2pairs = input_dirs.flatten().map {
                                        tuple(sample_names[it.toString()],
                                              file("${it}/hicpro/mapping/*bwt2pairs.bam", type: "file", checkIfExists: true))
                                     }
    input_validpairs = input_dirs.flatten().map {
                                        tuple(sample_names[it.toString()],
                                              file("${it}/hicpro/valid_pairs/*.allValidPairs", type: "file", checkIfExists: true))
                                      }
    input_stats = input_dirs.flatten().map {
                                    tuple(sample_names[it.toString()],
                                          file("${it}/hicpro/stats/", type: "dir", checkIfExists: true),
                                          file("${it}/hicpro/valid_pairs/stats/*.mergestat", type: "file", checkIfExists: true))
                                      }
    
    coolers_by_condition = input_coolers.map { tuple(it[0], it[2]) }
                                        .groupTuple(by: 0, size: 2)
                                        .map { tuple(it[0], it[1].flatten()) }
                      
    coolers_by_sample = input_coolers.map { tuple(it[1], it[2].flatten()) }

    multiqc(input_dirs)

    generate_blacklist(file(params.assembly_gaps),
                       file(params.cytoband))

    cooler_merge(coolers_by_condition)

    cooler_to_hic(coolers_by_sample.mix(cooler_merge.out.cool))

    cooler_zoomify(coolers_by_sample.mix(cooler_merge.out.cool),
                   generate_blacklist.out.bed)

    compress_bwt2pairs(input_bwt2pairs,
                       file(params.fasta))

    compress_validpairs(input_validpairs)

    compress_stats(input_stats)

    // NOTE label and sample mappings is guaranteed to be correct
    // regardless of the order in which items are emitted by compress_stats()
    // because plot_stats is globbing file names, thus sorting files by their
    // sample number
    plot_stats(compress_stats.out.tar.map { it[1] }.collect(),
               params.plot_pretty_labels)
}

process multiqc {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: {
                                        label = file(folder).getName()
                                        "${label}/${it}"
                                        }
    input:
        path folder

    output:
        path "*.html", emit: html
        path "*.tar.gz", emit: tar

    shell:
        '''
        multiqc .

        tar -cf - multiqc_data/ | gzip -9 > multiqc_data.tar.gz
        '''
}

process generate_blacklist {
    input:
        path assembly_gaps
        path cytoband

    output:
        path "*.bed", emit: bed

    shell:
        '''
        set -o pipefail

        cat <(zcat '!{assembly_gaps}' | cut -f 2-) \
            <(zcat '!{cytoband}' | grep 'acen$') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 > blacklist.bed
        '''
}

process cooler_merge {
    input:
        tuple val(label), path(cool)
        
    output:
        tuple val("${label}_merged"), path("*.cool"), emit: cool
        
    shell:
        '''
        cooler merge '!{label}_merged.cool' *.cool
        '''
}

process cooler_to_hic {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${label}/${label}.hic" }

    label 'process_medium'
    label 'process_high_memory'

    input:
        tuple val(label), path(cool)

    output:
        path "*.hic", emit: hic

    shell:
        memory_mb=task.memory.toMega() - 750
        '''
        set -o pipefail
        trap 'rm -rf tmp/' EXIT

        mkdir tmp/

        zstd_wrapper='tmp/zstd_wrapper.sh'

        printf '%s\\n%s\\n' \
            '#!/usr/bin/env bash' \
            'if [[ $# == 0 ]]; then zstd --adapt -T!{task.cpus}; else zstd -d; fi' > "$zstd_wrapper"
        chmod 755 "$zstd_wrapper"

        cat "$zstd_wrapper"

        cooler dump -t chroms '!{cool}' > tmp/chrom.sizes

        # See here for the format produced by AWK command
        # https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
        # IMPORTANT: juicer_tools expect chromosome names to be sorted, so that's why we need the sort command
        # Unfortunately pipes are not supported, so we have to write pixels to a temporary file

        cooler dump -t pixels --join '!{cool}' |
            awk -F '\\t' 'BEGIN{ OFS=FS } {print "1",$1,($2+$3)/2,"0","0",$4,($5+$6)/2,"1",$7}' |
            sort -k2,2 -k6,6 -T tmp/ -S !{memory_mb}M --compress-program="$zstd_wrapper" --parallel=!{task.cpus} |
            pigz -9 -p !{task.cpus} > tmp/pixels.txt.gz

        java -Xms512m -Xmx!{memory_mb}m -jar \
            /usr/local/share/java/juicer_tools/*.jar \
            pre -j !{task.cpus} tmp/pixels.txt.gz matrix.hic tmp/chrom.sizes
        '''
}

process cooler_zoomify {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${label}/${label}.mcool" }                                             

    label 'process_medium'
    label 'process_long'

    input:
        tuple val(label), path(cool)
        path blacklist
    
    output:
        tuple val(label), path("*.mcool"), emit: mcool

    shell:
        '''
        balance_args='-p !{task.cpus} --blacklist '!{blacklist}' --max-iters 500'

        cooler zoomify -p !{task.cpus}                  \
                       -r N                             \
                       --balance                        \
                       --balance-args="${balance_args}" \
                       '!{cool}'
        '''
}

process compress_bwt2pairs {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${label}/${it}" }                                             
    
    label 'process_long'
    label 'process_high'

    input:
        tuple val(label), path(bams)
        path reference_fna

    output:
        tuple val(label), path("*.cram"), emit: cram

    shell:
        '''
        set -o pipefail

        samtools merge -c -r -u -@!{task.cpus} -o - *.bam |
        samtools sort -u -@!{task.cpus}                   |
        samtools view --reference '!{reference_fna}'      \
                      --output-fmt cram,archive,use_lzma  \
                      -@!{task.cpus}                      \
                      -o '!{label}.bwt2pairs.cram'
        '''
}


process compress_validpairs {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${label}/${it}" }                                             

    label 'process_long'
    label 'process_high'

    input:
        tuple val(label), path(validpairs)

    output:
        tuple val(label), path("*.xz"), emit: xz

    shell:
        '''
        xz -T!{task.cpus} -9 --extreme -k -f '!{validpairs}'
        '''
}

process compress_stats {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${label}/stats.tar.xz" }

    label 'process_short'

    input:
        tuple val(label), path(stats_dir), path(valid_pairs_stats)

    output:
        tuple val(label), path("*.xz"), emit: tar

    shell:
        outprefix="${label}_stats"
        '''
        set -o pipefail

        mkdir '!{outprefix}'

        rsync -aPLv '!{stats_dir}/'* '!{valid_pairs_stats}' '!{outprefix}/'

        tar -chf - '!{outprefix}' |
            xz -T!{task.cpus} -9 --extreme > '!{outprefix}.tar.xz'
        '''
}

process plot_stats {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path archives
        val pretty_labels

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        '''
        generate_hic_mapping_report.py \
            *.tar.xz \
            hic_mapping_report \
            --sample-labels '!{pretty_labels}'

        mkdir plots
        mv *.png *.svg plots/
        '''
}
