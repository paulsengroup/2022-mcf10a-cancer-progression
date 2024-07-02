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

    Channel.fromPath(params.domains_cis_merged, checkIfExists: true)
        .collect()
        .set { domains }

    pair_maximal_cliques(
        cliques,
        labels_merged
    )

    plot_maximal_clique_alluvials(
        pair_maximal_cliques.out.tsv
    )

    plot_median_clique_genomic_span(
        cliques,
        domains,
        labels_merged
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
            *{WT,T1,C1}_cis_cliques.tsv.gz  \\
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

process plot_median_clique_genomic_span {
    publishDir "${params.output_dir}/plots/cliques/", mode: 'copy'

    input:
        path cliques
        path domains
        val labels

    output:
        path "*.svg", emit: svg

    shell:
        cliques_str=cliques.join(" ")
        domains_str=domains.join(" ")
        '''
        plot_median_clique_genomic_span_distribution.py \\
            --cliques *{WT,T1,C1}_cis_cliques.tsv.gz  \\
            --domains *{WT,T1,C1}_cis_domains.bed.gz  \\
            --labels='!{labels}' \\
            -o cis_clique_span_distribution.svg
        '''
}
