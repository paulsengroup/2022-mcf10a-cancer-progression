#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


workflow {

    labels_merged = params.labels_merged.join(",")

    Channel.fromPath(params.cliques_cis_merged, checkIfExists: true)
        .collect()
        .set { cliques }

    pair_maximal_cliques(
        cliques,
        labels_merged
    )

    plot_maximal_clique_alluvials(
        pair_maximal_cliques.out.tsv
    )

}

process pair_maximal_cliques {
    input:
        path cliques
        val labels

    output:
        path "*.tsv", emit: tsv

    shell:
        cliques_str=cliques.join(" ")
        '''
        pair_maximal_cliques.py  \\
            !{cliques_str}       \\
            --labels='!{labels}' > paired_cliques.tsv
        '''
}

process plot_maximal_clique_alluvials {
    publishDir "${params.output_dir}/plots/cliques/", mode: 'copy'

    input:
        path paired_cliques

    output:
        path "*.svg", emit: svg

    shell:
        '''
        plot_clique_alluvials.py                                    \\
            '!{paired_cliques}'                                     \\
            --path-to-plotting="$(which plot_clique_alluvials.r)"   \\
            -o 'cis_alluvial'
        '''
}
