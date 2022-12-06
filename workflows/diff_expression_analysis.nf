#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    preprocess_count_matrix(Channel.of(
                                tuple("genes", file(params.gene_count_matrix)),
                                tuple("transcripts", file(params.transcript_count_matrix))),
                            file(params.control_to_sample_mappings))

    run_deseq2(preprocess_count_matrix.out)

    run_elixir_gost(run_deseq2.out.map { it[1] }
                              .flatten(),
                    params.elixir_gost_lfc_lb,
                    params.elixir_gost_lfc_ub,
                    params.elixir_gost_pval)
}


process preprocess_count_matrix {
    input:
        tuple val(label), path(gene_count_matrix)
        path sample_mappings

    output:
        tuple val(label), path("*.tsv"), env(CONTRAST)

    shell:
        '''
        '!{params.script_dir}/map_rnaseq_samples_to_controls.py' \
            --round                                              \
            '!{gene_count_matrix}'                               \
            '!{sample_mappings}'                                 \
            .

        CONTRAST="$(head -n 1 *.contrast | tr -d '\\n')"
        '''
}

process run_deseq2 {
    publishDir "${params.output_dir}/", mode: 'copy'

    input:
        tuple val(label), path(count_matrix), val(contrast)

    output:
        tuple val(label), path("*.tsv"), path("*.pdf")

    shell:
        '''
        '!{params.script_dir}/diff_expression_analysis.r' \
            --count_matrix='!{count_matrix}'              \
            --contrast=!{contrast}                        \
            -o '!{label}_'

        rm -f Rplots.pdf
        '''
}

process run_elixir_gost {
    publishDir "${params.output_dir}/", mode: 'copy'

    input:
        path de_genes
        val lfc_lb
        val lfc_ub
        val pval

    output:
        path "*.tsv"

    shell:
        outname="${de_genes.baseName}_gost.tsv"
        '''
        '!{params.script_dir}/run_elixir_gost.py' \
            --pval-cutoff '!{pval}'               \
            --lfc-cutoffs '!{lfc_lb}' '!{lfc_ub}' \
            < '!{de_genes}' > '!{outname}'
        '''
}
