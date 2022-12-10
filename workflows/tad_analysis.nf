#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    hicexplorer_find_tads(Channel.fromPath(params.mcools),
                          params.resolution)

    generate_tad_report(hicexplorer_find_tads.out.domains.collect(),
                        hicexplorer_find_tads.out.scores.collect(),
                        params.repl_pretty_labels,
                        params.cond_pretty_labels)
}

process hicexplorer_find_tads {
    publishDir "${params.output_dir}", mode: 'copy'

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

process generate_tad_report {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path domains
        path scores
        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        '''
        '!{params.script_dir}/generate_tad_report.py' \
            GRCh38_???_*_domains.bed.gz \
            --output-prefix=report_replicates \
            --labels='!{labels_replicates}'

        '!{params.script_dir}/generate_tad_report.py' \
            *{WT,T1,C1}_merged_domains.bed.gz \
            --output-prefix=report_conditions \
            --labels='!{labels_conditions}'

        '!{params.script_dir}/generate_insulation_report.py' \
            GRCh38_???_*_score.bedgraph.gz \
            --output-prefix=report_replicates_insulation \
            --labels='!{labels_replicates}'

        '!{params.script_dir}/generate_insulation_report.py' \
            *{WT,T1,C1}_merged_score.bedgraph.gz \
            --output-prefix=report_conditions_insulation \
            --labels='!{labels_conditions}'

        mkdir plots
        mv *.svg *.png plots/
        '''
}
