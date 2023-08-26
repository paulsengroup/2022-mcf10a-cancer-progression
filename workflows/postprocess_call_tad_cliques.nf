#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


workflow {
    Channel.empty()
        .mix(Channel.fromPath(params.cliques_cis, checkIfExists: true).collect().map { tuple("cis-only", it) })
        .mix(Channel.fromPath(params.cliques_trans, checkIfExists: true).collect().map { tuple("trans-only", it) })
        .mix(Channel.fromPath(params.cliques_all, checkIfExists: true).collect().map { tuple("all", it) })
        .map{ tuple(it[0], it[1], params.labels.join(","), "") }
        .set { cliques }

    Channel.empty()
        .mix(Channel.fromPath(params.cliques_cis_wt, checkIfExists: true).collect().map { tuple("cis-only", it) })
        .mix(Channel.fromPath(params.cliques_trans_wt, checkIfExists: true).collect().map { tuple("trans-only", it) })
        .mix(Channel.fromPath(params.cliques_all_wt, checkIfExists: true).collect().map { tuple("all", it) })
        .map{ tuple(it[0], it[1], params.labels.join(","), "aggregate_on_wt") }
        .set { cliques_wt }

    plot_maximal_clique_sizes(
        cliques.mix(cliques_wt)
    )
}

process plot_maximal_clique_sizes {
    publishDir "${params.output_dir}/${output_folder}/plots", mode: 'copy'


    tag "${type}"

    input:
        tuple val(type),
              path(cliques),
              val(labels),
              val(output_folder)

    output:
        tuple val(type),
              path("*.png"), emit: png

        tuple val(type),
              path("*.svg"), emit: svg

    shell:
        '''
        plot_maximal_clique_sizes.py            \\
            *{WT,T1,C1}_*cliques.tsv.gz         \\
            -o '!{type}_maximal_clique_size_abs'

        plot_maximal_clique_sizes.py            \\
            *{WT,T1,C1}_*cliques.tsv.gz         \\
            --stat='density'                    \\
            -o '!{type}_maximal_clique_size_rel'
        '''
}

