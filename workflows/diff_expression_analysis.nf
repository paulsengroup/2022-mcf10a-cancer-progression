#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    run_deseq2(file(params.count_matrix, checkIfExists: true),
               file(params.design_table, checkIfExists: true),
               Channel.of(params.lfc_cutoffs).flatten(),
               params.fdr_alpha)
}


process run_deseq2 {
    publishDir "${params.output_dir}/", mode: 'copy',
                                        saveAs: { "lfc_${lfc_cutoff}/de_genes.tsv.gz" }

    input:
        path count_matrix
        path design_table
        val lfc_cutoff
        val fdr_alpha

    output:
        path "*.tsv.gz"

    shell:
        '''
        set -o pipefail

        diff_expression_analysis.py      \\
            '!{count_matrix}'            \\
            '!{design_table}'            \\
            --lfc-thresh='!{lfc_cutoff}' \\
            --fdr-alpha='!{fdr_alpha}' | gzip -9 > de_genes.tsv.gz
        '''
}
