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
    
    cooler_zoomify(coolers_by_sample.mix(cooler_merge.out.cool),
                   generate_blacklist.out.bed)

    compress_bwt2pairs(input_bwt2pairs,
                       file(params.fasta))

    compress_validpairs(input_validpairs)

    compress_stats(input_stats)
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
                                       saveAs: { "${label}/${it}" }

    label 'process_short'

    input:
        tuple val(label), path(stats_dir), path(valid_pairs_stats)

    output:
        tuple val(label), path("*.xz"), emit: xz

    shell:
        outprefix="${label}_stats"
        '''
        set -o pipefail

        mkdir '!{outprefix}'

        rsync -aPLv '!{stats_dir}/'* '!{valid_pairs_stats}' '!{outprefix}/'

        tar -chf - '!{outprefix}' |
            xz -T!{task.cpus} -9 --extreme > 'stats.tar.xz'
        '''
}