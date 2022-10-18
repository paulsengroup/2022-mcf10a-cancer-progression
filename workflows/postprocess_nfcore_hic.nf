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
                                          file("${it}/hicpro/stats/", type: "dir", checkIfExists: true))
                                 }
    
    coolers_by_condition = input_coolers.map { tuple(it[0], it[2]) }
                                        .groupTuple(by: 0, size: 2)
                                        .map { tuple(it[0], it[1].flatten()) }
                      
    coolers_by_sample = input_coolers.map { tuple(it[1], it[2].flatten()) }


    cooler_merge(coolers_by_condition)
    
    cooler_zoomify(coolers_by_sample.mix(cooler_merge.out.cool))

    compress_bwt2pairs(input_bwt2pairs,
                       file(params.fasta))

    compress_validpairs(input_validpairs)

    compress_stats(input_stats)
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
    
    output:
        tuple val(label), path("*.mcool"), emit: mcool

    shell:
        '''
        cooler zoomify -p !{task.cpus}                  \
                       -r N                             \
                       --balance                        \
                       --balance-args='-p !{task.cpus}' \
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
        tuple val(label), path(stats_dir)

    output:
        tuple val(label), path("*.xz"), emit: xz

    shell:
        '''
        set -o pipefail

        tar -cf - '!{stats_dir}' |
            xz -T!{task.cpus} -9 --extreme > 'stats.tar.xz'
        '''
}