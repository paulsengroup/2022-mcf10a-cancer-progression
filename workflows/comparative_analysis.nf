#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.of(params.resolutions)
       .flatten()
       .map { res ->
              def template = params.subcompartment_bedgraph_template as String
              def bgraph = template.replace('{{resolution}}', "${res}")

              tuple(res, file(bgraph, checkIfExists: true, glob: false, type: 'file'))
            }
       .set { subcomp_bgs }

    Channel.fromPath(["${params.nfcore_chip_folder}/*",
                      "${params.lads_folder}"],
            type: "dir",
            checkIfExists: true)
           .set { epigen_marker_folders }

    coverage_only_flags = [true, false]

    overlap_subcompartments_with_epigenetic_markers(
        subcomp_bgs.combine(coverage_only_flags),
        epigen_marker_folders.collect(),
        file(params.subcomp_marker_file_table, checkIfExists: true),
        params.condition_labels?.join(',')
    )

    plot_subcompartment_vs_epigenetic_markers(
        overlap_subcompartments_with_epigenetic_markers.out.pickle
    )

    annotate_domains_with_subcompartments(
        subcomp_bgs,
        file(params.tads, checkIfExists: true),
        file(params.cliques, checkIfExists: true),
        params.subcomp_clique_labels
    )

    annotate_domains_with_subcompartments.out.tads
        .groupTuple()
        .map { tuple(it[0], 'tad', it.drop(1).flatten()) }
        .set { annotated_tads }

    annotate_domains_with_subcompartments.out.cliques
        .groupTuple()
        .map { tuple(it[0], 'clique', it.drop(1).flatten()) }
        .set { annotated_cliques }

    cluster_domains_by_subcompartment(
        annotated_tads.mix(annotated_cliques),
        params.repl_labels.join(','),
        params.hdbscan_min_cluster_size,
        params.hdbscan_min_samples,
        params.hdbscan_cluster_selection_method,
        params.hdbscan_dist
    )

    cluster_domains_by_subcompartment.out.pickle
        .join(cluster_domains_by_subcompartment.out.tsv,
              by: [0, 1],
              failOnDuplicate: true,
              failOnMismatch: true)
        .combine(params.repl_labels + ['all'])
        .set { plot_domain_clusters_inputs }
        // [resolution, label, clusterer, cluster, condition]

    plot_domain_clusters(
        plot_domain_clusters_inputs
    )

    plot_clique_subcomp_composition(
        annotated_cliques,
        params.repl_labels.join(',')
    )

    overlap_subcomps_with_expression_lvls(
        subcomp_bgs,
        file(params.expression_table_tpm, checkIfExists: true),
        file(params.gtf, checkIfExists: true),
        file(params.rnaseq_design_table, checkIfExists: true),
        params.condition_labels.join(',')
    )

    overlap_subcomps_with_de_genes_heatmaps(
        subcomp_bgs,
        file(params.deg_table, checkIfExists: true),
        file(params.gtf, checkIfExists: true)
    )
}

process overlap_subcompartments_with_epigenetic_markers {
    publishDir "${params.output_dir}/subcomps_vs_epigenetic_markers/${resolution}", mode: 'copy'
    label 'process_medium'

    input:
        tuple val(resolution),
              path(subcompartments),
              val(coverage_only)
        path nfcore_chip_folder
        path marker_file_table
        val pretty_labels

    output:
        tuple val(resolution), path("*.pickle"), emit: pickle
        tuple val(resolution), path("*.tsv.gz"), emit: tsv

    shell:
        outprefix = subcompartments.getSimpleName() + (coverage_only ? "_coverage" : "")
        coverage_only_flag = coverage_only ? "--coverage-only" : ""
        '''
        overlap_subcomps_with_epigen_markers.py \\
            '!{subcompartments}' \\
            '!{marker_file_table}' \\
            '!{outprefix}' \\
            --labels '!{pretty_labels}' \\
            --nproc=!{task.cpus} \\
            !{coverage_only_flag}

        pigz -9 -p !{task.cpus} '!{outprefix}.tsv'
        '''
}

process plot_subcompartment_vs_epigenetic_markers {
    publishDir "${params.output_dir}/subcomps_vs_epigenetic_markers/", mode: 'copy'
    label 'very_short'

    input:
        tuple val(resolution),
              path(pickled_df)

    output:
        tuple val(resolution), path("${resolution}/plots/*.png"), emit: png
        tuple val(resolution), path("${resolution}/plots/*.svg"), emit: svg

    shell:
        outprefix="${resolution}/plots/${pickled_df.simpleName}"
        '''
        plot_subcomp_epigen_marker_overlaps.py '!{pickled_df}' '!{outprefix}'
        '''
}


process annotate_domains_with_subcompartments {
    publishDir "${params.output_dir}/annotated_domains/", mode: 'copy'
    label 'process_short'

    input:
        tuple val(resolution),
              path(subcompartments)

        path tads
        path cliques
        val domain_names

    output:
        tuple val(resolution), path("${resolution}/*.bed.gz"), emit: tads
        tuple val(resolution), path("${resolution}/*.tsv.gz"), emit: cliques

    shell:
        outprefix="${resolution}/${subcompartments.simpleName}"
        domain_names_str=domain_names.join(" ")
        '''
        annotate_domains_with_subcompartments.py \\
            '!{subcompartments}' \\
            --domains *MCF10A_{WT,T1,C1}_cis_domains.bed.gz \\
            --cliques *MCF10A_{WT,T1,C1}_cis_cliques.tsv.gz \\
            --domain-names !{domain_names_str} \\
            --output-folder '!{resolution}'
        '''
}

process cluster_domains_by_subcompartment {
    publishDir "${params.output_dir}/annotated_domains/", mode: 'copy',
                                                          saveAs: { "${resolution}/clusters/${it}" }
    label 'process_short'

    input:
        tuple val(resolution),
              val(label),
              path(annotated_domains)

        val labels
        val min_size
        val min_samples
        val selection_method
        val metric

    output:
        tuple val(resolution),
              val(label),
              path("*.pickle.*"),
        emit: pickle
        tuple val(resolution),
              val(label),
              path("*.tsv.gz"),
        emit: tsv

    shell:
        outprefix = label
        '''
        cluster_domains_by_subcompartment_state.py \\
            --labels '!{labels}' \\
            --min-cluster-size '!{min_size}' \\
            --min-samples '!{min_samples}' \\
            --cluster-selection-method '!{selection_method}' \\
            --distance-metric '!{metric}' \\
            *MCF10A_{WT,T1,C1}*.gz \\
            --output-prefix='!{outprefix}'
        '''
}

process plot_domain_clusters {
    publishDir "${params.output_dir}/annotated_domains/", mode: 'copy',
                                                          saveAs: { "${resolution}/clusters/plots/${it}" }
    label 'process_short'

    input:
        tuple val(resolution),
              val(label),
              path(clusterer),
              path(clusters),
              val(condition_of_interest)

    output:
        tuple val(resolution),
              val(label),
              val(condition_of_interest),
              path("*.png"),
        emit: png
        tuple val(resolution),
              val(label),
              val(condition_of_interest),
              path("*.svg"),
        emit: svg

    shell:
        '''
        plot_domain_subcompartment_clusters.py \\
            '!{clusters}' \\
            '!{clusterer}' \\
            --plot-type=line \\
            --plot-title='!{label} ('!{condition_of_interest}')' \\
            --label-to-highlight='!{condition_of_interest}' \\
            --output-prefix='!{label}_!{condition_of_interest}'

        plot_domain_subcompartment_clusters.py \\
            '!{clusters}' \\
            '!{clusterer}' \\
            --plot-type=barplot \\
            --plot-title='!{label} ('!{condition_of_interest}')' \\
            --label-to-highlight='!{condition_of_interest}' \\
            --output-prefix='!{label}_!{condition_of_interest}'
        '''
}

process plot_clique_subcomp_composition {
    publishDir "${params.output_dir}/annotated_domains/", mode: 'copy',
                                                          saveAs: { "${resolution}/clusters/plots/${it}" }
    label 'process_short'

    input:
        tuple val(resolution),
              val(label),
              path(annotated_domains)

        val labels

    output:
        tuple val(resolution),
              val(label),
              path("*.svg"),
        emit: svg
        tuple val(resolution),
              val(label),
              path("*.png"),
        emit: png

    shell:
        outprefix = label
        '''
        plot_clique_subcomp_composition.py \\
            *{WT,T1,C1}*.tsv.gz \\
            --labels !{labels} \\
            --output-prefix='!{outprefix}'
        '''
}

process overlap_subcomps_with_expression_lvls {
    publishDir "${params.output_dir}/subcomps_vs_rnaseq/plots/", mode: 'copy'
    label 'process_short'

    input:
        tuple val(resolution),
              path(subcomps)
        path expression_table
        path gtf
        path sample_name_mappings
        val labels

    output:
        path "${resolution}/*.png", emit: png
        path "${resolution}/*.svg", emit: svg

    shell:
        '''
        mkdir '!{resolution}'
        overlap_subcomps_with_expression.py \\
            '!{subcomps}' \\
            '!{expression_table}' \\
            '!{gtf}' \\
            '!{resolution}/subcomps_vs_expression' \\
            --labels='!{labels}' \\
            --sample-name-mappings-tsv=<(tail -n +2 '!{sample_name_mappings}')
            # tail -n +2 is used to skip the header
        '''
}

process overlap_subcomps_with_de_genes_heatmaps {
    publishDir "${params.output_dir}/subcomps_vs_rnaseq/plots/", mode: 'copy'
    label 'process_short'

    input:
        tuple val(resolution),
              path(subcomps)
        path deg_table
        path gtf

    output:
        path "${resolution}/*.png", emit: png
        path "${resolution}/*.svg", emit: svg

    shell:
        '''
        mkdir '!{resolution}'

        overlap_subcomps_with_de_genes.py \\
            '!{subcomps}' \\
            MCF10A_WT_vs_MCF10A_T1.tsv.gz \\
            '!{gtf}' \\
            '!{resolution}/MCF10A_WT_vs_MCF10A_T1_subcomps_vs_deg' \\
            --contrast MCF10A_WT \\
            --condition MCF10A_T1 \\
            --plot-type=heatmap \\
            --padj 0.01 \\
            --lfc 2

        overlap_subcomps_with_de_genes.py \\
            '!{subcomps}' \\
            MCF10A_WT_vs_MCF10A_C1.tsv.gz \\
            '!{gtf}' \\
            '!{resolution}/MCF10A_WT_vs_MCF10A_C1_subcomps_vs_deg' \\
            --contrast MCF10A_WT \\
            --condition MCF10A_C1 \\
            --plot-type=heatmap \\
            --padj 0.01 \\
            --lfc 2

        overlap_subcomps_with_de_genes.py \\
            '!{subcomps}' \\
            MCF10A_T1_vs_MCF10A_C1.tsv.gz \\
            '!{gtf}' \\
            '!{resolution}/MCF10A_T1_vs_MCF10A_C1_subcomps_vs_deg' \\
            --contrast MCF10A_T1 \\
            --condition MCF10A_C1 \\
            --plot-type=heatmap \\
            --padj 0.01 \\
            --lfc 2
        '''
}
