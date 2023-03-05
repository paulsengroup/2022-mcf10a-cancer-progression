#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {

    Channel.fromPath(params.subcompartment_bedgraph, checkIfExists: true)
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
        params.condition_labels
    )

    plot_subcompartment_vs_epigenetic_markers(
        overlap_subcompartments_with_epigenetic_markers.out.pickle
    )
}

process overlap_subcompartments_with_epigenetic_markers {
    publishDir "${params.output_dir}/subcomps_vs_epigenetic_markers/", mode: 'copy'
    label 'process_medium'

    input:
        tuple path(subcompartments),
              val(coverage_only)
        path nfcore_chip_folder
        path marker_file_table
        val pretty_labels

    output:
        path "*.pickle", emit: pickle
        path "*.tsv.gz", emit: tsv

    shell:
        outprefix = subcompartments.getSimpleName() + (coverage_only ? "_coverage" : "")
        coverage_only_flag = coverage_only ? "--coverage-only" : ""
        '''
        ls -lah .
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
    publishDir "${params.output_dir}/subcomps_vs_epigenetic_markers/plots/", mode: 'copy'
    label 'very_short'

    input:
        path pickled_df

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg

    shell:
        outprefix = pickled_df.getSimpleName()
        '''
        plot_subcomp_epigen_marker_overlaps.py '!{pickled_df}' '!{outprefix}'
        '''
}
