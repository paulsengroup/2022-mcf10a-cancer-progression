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

    run_go_figure(run_elixir_gost.out.tsv.flatten())
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
        tuple val(label), path("*.tsv"), path("plots/*.pdf")

    shell:
        '''
        '!{params.script_dir}/diff_expression_analysis.r' \
            --count_matrix='!{count_matrix}'              \
            --contrast=!{contrast}                        \
            -o '!{label}_'

        rm -f Rplots.pdf

        mkdir plots
        mv *.pdf plots/
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
        path "*.tsv", emit: tsv

    shell:
        outprefix="${de_genes.baseName}_gost"
        '''
        lfc_lb=($(echo '!{lfc_lb}' | tr ',' ' '))
        lfc_ub=($(echo '!{lfc_ub}' | tr ',' ' '))
        pval=($(echo '!{pval}' | tr ',' ' '))

        for i in "${!lfc_lb[@]}"; do
            pv="${pval[$i]}"
            lb="${lfc_lb[$i]}"
            ub="${lfc_ub[$i]}"

            '!{params.script_dir}/run_elixir_gost.py' \
                --pval-cutoff "$pv"                   \
                --lfc-cutoffs "$lb" "$ub"             \
                < '!{de_genes}' > "!{outprefix}_${lb}_${ub}_${pv}.tsv"
        done
        '''
}

process run_go_figure {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    input:
        path functional_annotation

    output:
        path "*.png", emit: png

    shell:
        outprefix="${functional_annotation.baseName}"
        '''
        set -o pipefail

        trap 'rm -f annotation.tsv' EXIT
        grep '^GO:' '!{functional_annotation}' | cut -f 2,4 > annotation.tsv

        gofigure.py --input annotation.tsv --font_size=xx-small --description_limit=50 --max_label=20 --output '!{outprefix}'

        mv '!{outprefix}/'* .
        '''
}