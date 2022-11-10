#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    preprocess_count_matrix(file(params.gene_count_matrix),
                            file(params.control_to_sample_mappings))

    run_deseq2(preprocess_count_matrix.out.tsv.flatten(),
               preprocess_count_matrix.out.contrast.flatten())
}


process preprocess_count_matrix {
    input:
        path gene_count_matrix
        path sample_mappings

    output:
        path "*.tsv", emit: tsv
        path "*.contrast", emit: contrast

    shell:
        '''
        '!{params.script_dir}/map_rnaseq_samples_to_controls.py' \
            --round                                              \
            '!{gene_count_matrix}'                               \
            '!{sample_mappings}'                                 \
            .
        '''
}

process run_deseq2 {
    publishDir "${params.output_dir}/", mode: 'copy'

    input:
        path gene_count_matrix
        path contrast

    output:
        path "*.tsv", emit: tsv
        path "*.pdf", emit: pdf
        env CONTRAST, emit: contrast

    shell:
        '''
        CONTRAST="$(head -n 1 '!{contrast}' | tr -d '\\n')"

        '!{params.script_dir}/diff_expression_analysis.r' \
            --count_matrix='!{gene_count_matrix}'         \
            --contrast="$CONTRAST"                        \
            -o .

        rm -f Rplots.pdf
        '''
}
