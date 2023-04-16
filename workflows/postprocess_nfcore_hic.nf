#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


def collect_files(prefix, sample_id, suffix, type = "file") {
    def files = file("${prefix}/${sample_id}*${suffix}",
                     type: type,
                     checkIfExists: true)

    def condition_id = sample_id.replaceAll(/_REP\d+$/, "")

    return tuple(sample_id, condition_id, files)
}


workflow {
    sample_sheet.splitCsv(sep: ",", header: true)
        .map { it.sample }
        .unique()
        .set { sample_ids }

    input_dir = file(params.input_dir, type: "dir", checkIfExists: true)

    sample_ids
        .map { collect_files("${input_dir}/hicpro/mapping", it, ".bwt2pairs.bam") }
        .set { input_bwt2pairs }

    sample_ids
        .map { collect_files("${input_dir}/hicpro/valid_pairs", it, ".allValidPairs") }
        .set { input_valid_pairs }

    sample_ids
        .map { collect_files("${input_dir}/hicpro/stats", it, "/", "dir") }
        .set { input_stats }

    cooler_cload(input_validpairs,
                 file(params.chrom_sizes),
                 params.assembly_name,
                 params.cooler_base_resolution)

    cooler_cload.out.cool.set { coolers_by_sample }

    // Group coolers by condition
    coolers_by_sample
        .groupTuple(by: 1, size: 2)
        // [condition_id, [sample_ids], [coolers]]
        .set { coolers_by_condition }

    // multiqc(input_dirs)  TODO: update me!

    generate_blacklist(file(params.assembly_gaps),
                       file(params.cytoband))

    cooler_merge(coolers_by_condition)

    Channel.empty()
        .mix(coolers_by_sample.map { tuple(it[0], it[2]) },     // [sample_id, cooler]
             coolers_by_condition.map { tuple(it[0], it[2] ) }) // [condition_id, cooler]
        .set { coolers }

    cooler_to_hic(coolers)

    cooler_zoomify(coolers)
    cooler_balance(cooler_zoomify.out.mcool,
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


process archive_multiqc {
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

        cat <(zcat '!{assembly_gaps}' | cut -f 2-) \\
            <(zcat '!{cytoband}' | grep 'acen$') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 > blacklist.bed
        '''
}


process cooler_cload {
    input:
        tuple val(sample),
              val(condition),
              path(pairs)

        val chrom_sizes
        val assembly
        val resolution

    output:
        tuple val(sample),
              val(condition),
              path("*.cool"),
        emit: cool

    shell:
        outname = "${pairs.simpleName}.cool"
        '''
        chrom_sizes="$(mktemp chrom.sizes.XXXXXX)"

        grep '^chr[XY0-9]\\+[[:space:]]' '!{chrom_sizes}' > "$chrom_sizes"

        cooler cload \\
            "$chrom_sizes:!{resolution}" \\
            '!{pairs}'                   \\
            '!{outname}'                 \\
            --assembly='!{assembly}'
        '''
}

process cooler_merge {
    input:
        tuple val(condition),
              val(samples),
              path(coolers)

    output:
        tuple val("${condition}_merged"),
              val(samples),
              path("*.cool"),
        emit: cool

    shell:
        outname = "${condition}_merged.cool"
        '''
        cooler merge '!{outname}' *.cool
        '''
}

process cooler_to_hic {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "hic/${label}.hic" }

    label 'process_medium'

    input:
        tuple val(abel),
              path(cool)

    output:
        tuple val(label),
              path("*.hic"),
        emit: hic

    shell:
        memory_gb=task.memory.toGiga()
        out="${cool.baseName}.hic"
        '''
        trap 'rm -rf tmp/' EXIT

        hic_tools_jar=(/usr/local/share/java/hic_tools/hic_tools*.jar)
        juicer_tools_jar=(/usr/local/share/java/juicer_tools/juicer_tools*.jar)

        cool2hic-ng \\
            --hic-tools-jar "${hic_tools_jar[*]}"       \\
            --juicer-tools-jar "${juicer_tools_jar[*]}" \\
            --tmpdir=tmp/                               \\
            '!{cool}'                                   \\
            '!{out}'                                    \\
            --nproc '!{task.cpus}'                      \\
            -Xmx '!{memory_gb}g'                        \\
            --resolutions 1000 2000 5000 10000 20000    \\
                          50000 100000 200000 500000    \\
                          1000000 2000000 5000000 10000000
        '''
}

process cooler_zoomify {
    label 'process_medium'
    label 'process_long'

    input:
        tuple val(label),
              path(cool)

    output:
        tuple val(label),
              path("*.mcool"),
        emit: mcool

    shell:
        '''
        cooler zoomify -p !{task.cpus} \\
                       -r N            \\
                       '!{cool}'
        '''
}

process cooler_balance {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "mcools/${label}.mcool" }

    label 'process_medium'
    label 'process_long'

    input:
        tuple val(label),
              path(mcool)
        path blacklist

    output:
        tuple val(label),
              path("*.mcool.new"),
        emit: mcool

    shell:
        memory_gb=task.memory.toGiga()
        '''
        hic_tools_jar=(/usr/local/share/java/hic_tools/hic_tools*.jar)
        juicer_tools_jar=(/usr/local/share/java/juicer_tools/juicer_tools*.jar)

        cp '!{mcool}' '!{mcool}.new'
        cooler_balance.py \\
            --hic-tools-jar "${hic_tools_jar[*]}"       \\
            --juicer-tools-jar "${juicer_tools_jar[*]}" \\
            --nproc !{task.cpus}                        \\
            --blacklist '!{blacklist}'                  \\
            -Xmx !{memory_gb}g                          \\
            '!{mcool}'                                  \\
            '!{mcool}.new'
        '''
}

process compress_bwt2pairs {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "alignments/${it}" }

    label 'process_long'
    label 'process_high'

    input:
        tuple val(label),
              path(bams)
        path reference_fna

    output:
        tuple val(label),
              path("*.cram"),
        emit: cram

    shell:
        '''
        set -o pipefail

        samtools merge -c -r -u -@!{task.cpus} -o - *.bam |
        samtools sort -u -@!{task.cpus}                   |
        samtools view --reference '!{reference_fna}'      \\
                      --output-fmt cram,archive,use_lzma  \\
                      -@!{task.cpus}                      \\
                      -o '!{label}.bwt2pairs.cram'
        '''
}

process compress_validpairs {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "pairs/${it}" }

    label 'process_long'
    label 'process_high'

    input:
        tuple val(label),
              path(validpairs)

    output:
        tuple val(label),
              path("*.xz"),
        emit: xz

    shell:
        '''
        xz -T!{task.cpus} -9 --extreme -k -f '!{validpairs}'
        '''
}

process compress_stats {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "stats/stats.tar.xz" }

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
