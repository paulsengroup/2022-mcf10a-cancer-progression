#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    cooltools_insulation(Channel.fromPath(params.mcools),
                         params.insulation_windows,
                         params.resolution)

    hicexplorer_find_tads(Channel.fromPath(params.mcools),
                          params.resolution)
}

process cooltools_insulation {
    publishDir "${params.output_dir}/insulation", mode: 'copy'

    input:
        path cooler
        val windows
        val resolution

    output:
        path "*.tsv.gz", emit: tsv
        path "*.bw", emit: bw

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        cooltools insulation --bigwig               \
                             "$cooler"              \
                             !{windows}             \
                             -o '!{outprefix}.tsv'

        gzip -9 '!{outprefix}.tsv'
        '''
}

process hicexplorer_find_tads {
    publishDir "${params.output_dir}/tads", mode: 'copy'

    label 'process_medium'

    input:
        path cooler
        val resolution

    output:
        path "*boundaries.bed.gz", emit: boundaries
        path "*domains.bed.gz", emit: domains
        path "*score.bedgraph.gz", emit: scores
        path "*tad_score.bm.gz", emit: tad_scores

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        hicFindTADs -p '!{task.cpus}'          \
                    --matrix "$cooler"         \
                    --outPrefix '!{outprefix}' \
                    --correctForMultipleTesting fdr

        # Compress text files
        printf '%s\n' *.bed *.bedgraph *.bm |
            xargs -L 1 -P '!{task.cpus}' sh -c 'gzip -9 "$1"' sh
        '''
}