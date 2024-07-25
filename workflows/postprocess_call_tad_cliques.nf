#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


workflow {

    labels_merged = params.labels_merged.join(",")

    Channel.fromList(params.cliques_cis_merged)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { cliques_merged }

    Channel.fromList(params.domains_cis_merged)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { domains_merged }

    labels_repl = params.labels_repl.join(",")

    Channel.fromList(params.cliques_cis_repl)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { cliques_repl }

    Channel.fromList(params.domains_cis_repl)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { domains_repl }

    cliques_merged.join(domains_merged)
        .set { masking_tasks_merged }

    cliques_repl.join(domains_repl)
        .set { masking_tasks_repl }

    mask_cliques_merged(
        masking_tasks_merged,
        file(params.mask, checkIfExists: true)
    )

    mask_cliques_repl(
        masking_tasks_repl,
        file(params.mask, checkIfExists: true)
    )

    mask_cliques_merged.out.cliques
        .set { masked_cliques_merged }

    mask_cliques_repl.out.cliques
        .set { masked_cliques_repl }

    plot_maximal_tad_clique_size_merged(
        masked_cliques_merged.map { it[1] }.collect(),
        labels_merged
    )

    plot_maximal_tad_clique_size_repl(
        masked_cliques_repl.map { it[1] }.collect(),
        labels_repl
    )

    pair_maximal_cliques(
        cliques_merged.map{ it[1] }.collect(),
        labels_merged
    )

    plot_maximal_clique_alluvials(
        pair_maximal_cliques.out.tsv
    )

    plot_median_clique_genomic_span_merged(
        cliques_merged.map{ it[1] }.collect(),
        domains_merged.map{ it[1] }.collect(),
        labels_merged
    )

    plot_median_clique_genomic_span_repl(
        cliques_repl.map{ it[1] }.collect(),
        domains_repl.map{ it[1] }.collect(),
        labels_repl
    )
}

process mask_cliques_merged {
    input:
        tuple val(sample),
              path(cliques),
              path(domains)

        path mask

    output:
        tuple val(sample),
              path("*.tsv.gz"),
        emit: cliques

    shell:
    outname="${cliques.simpleName}.masked.tsv.gz"
    '''
    set -o pipefail

    mask_cliques.py '!{cliques}' '!{domains}' --mask='!{mask}' |
        gzip -9 > '!{outname}'
    '''
}

process mask_cliques_repl {
    input:
        tuple val(sample),
              path(cliques),
              path(domains)

        path mask

    output:
        tuple val(sample),
              path("*.tsv.gz"),
        emit: cliques

    shell:
    outname="${cliques.simpleName}.masked.tsv.gz"
    '''
    set -o pipefail

    mask_cliques.py '!{cliques}' '!{domains}' --mask='!{mask}' |
        gzip -9 > '!{outname}'
    '''
}

process plot_maximal_tad_clique_size_merged {
    publishDir "${params.output_dir_merged}/plots/cliques/", mode: 'copy'

    input:
        path cliques

        val labels

    output:
        path "*.svg", emit: svg

    shell:
        labels_str=labels.split(",").join(" ")
        '''
        plot_tad_maximal_clique_size_distribution.py \\
            *{WT,T1,C1}_cis_cliques.masked.tsv.gz \\
            --labels !{labels_str} \\
            -o 'cis_tad_max_clique_size_distribution_abs_with_masking'

        plot_tad_maximal_clique_size_distribution.py \\
            *{WT,T1,C1}_cis_cliques.masked.tsv.gz \\
            --labels !{labels_str} \\
            --stat='density' \\
            -o 'cis_tad_max_clique_size_distribution_rel_with_masking'
        '''
}

process plot_maximal_tad_clique_size_repl {
    publishDir "${params.output_dir_repl}/plots/cliques/", mode: 'copy'

    input:
        path cliques

        val labels

    output:
        path "*.svg", emit: svg

    shell:
        labels_str=labels.split(",").join(" ")
        '''
        plot_tad_maximal_clique_size_distribution.py \\
            *{WT,T1,C1}_{REP1,REP2}_cis_cliques.masked.tsv.gz \\
            --labels !{labels_str} \\
            -o 'cis_tad_max_clique_size_distribution_abs_with_masking'

        plot_tad_maximal_clique_size_distribution.py \\
            *{WT,T1,C1}_{REP1,REP2}_cis_cliques.masked.tsv.gz \\
            --labels !{labels_str} \\
            --stat='density' \\
            -o 'cis_tad_max_clique_size_distribution_rel_with_masking'
        '''
}

process pair_maximal_cliques {
    input:
        path cliques
        val labels

    output:
        path "*.tsv.gz", emit: tsv

    shell:
        '''
        pair_maximal_cliques.py  \\
            *{WT,T1,C1}_cis_cliques.tsv.gz  \\
            --labels='!{labels}' |
            gzip -9 > paired_cliques.tsv.gz
        '''
}

process plot_maximal_clique_alluvials {
    publishDir "${params.output_dir_merged}/plots/cliques/", mode: 'copy'

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

process plot_median_clique_genomic_span_merged {
    publishDir "${params.output_dir_merged}/plots/cliques/", mode: 'copy'

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

process plot_median_clique_genomic_span_repl {
    publishDir "${params.output_dir_repl}/plots/cliques/", mode: 'copy'

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
            --cliques *{WT,T1,C1}_{REP1,REP2}_cis_cliques.tsv.gz  \\
            --domains *{WT,T1,C1}_{REP1,REP2}_cis_domains.bed.gz  \\
            --labels='!{labels}' \\
            -o cis_clique_span_distribution.svg
        '''
}
