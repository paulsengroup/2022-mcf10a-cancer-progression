#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    run_deseq2(file(params.count_matrix, checkIfExists: true),
               file(params.design_table, checkIfExists: true),
               Channel.of(params.lfc_cutoffs).flatten())
}


process run_deseq2 {
    publishDir "${params.output_dir}/", mode: 'copy'

    input:
        path count_matrix
        path design_table
        val lfc_cutoff

    output:
        path "lfc_${lfc_cutoff}/*.tsv.gz", emit: tsv
        path "lfc_${lfc_cutoff}/*.rds", emit: res
        path "lfc_${lfc_cutoff}/dds.rds", emit: dds
        path "lfc_${lfc_cutoff}/r_sessioninfo.txt", emit: sessioninfo

    shell:
        '''
        run_deseq2.py           \\
            '!{count_matrix}'   \\
            '!{design_table}'   \\
            'lfc_!{lfc_cutoff}' \\
            --lfc-thresh='!{lfc_cutoff}'
        '''
}
