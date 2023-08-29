#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.mcools)
           .map { tuple(file(it).getSimpleName(), it) }
           .combine(params.resolutions)
            // label, resolution, cooler
           .map { tuple(it[0], it[2], it[1] )}
           .set{ coolers }


    compute_expected_cis(
        coolers
    )

    call_loops(
        compute_expected_cis.out.tsv
            .join(coolers, by: [0, 1]) // label, resolution
    )

    lowest_resolution = params.resolutions.max()
    merge_loops(
        call_loops.out.bedpe.groupTuple(size: params.resolutions.size()),
        lowest_resolution
    )

    detect_differential_loops(
        merge_loops.out.bedpe
            .map { it[1] }
            .collect(),
        lowest_resolution
    )

    generate_loop_report(
        Channel.fromPath(params.mcools)
            .collect(),
        merge_loops.out.bedpe
            .map { it[1] }
            .collect(),
        lowest_resolution
    )

}

process compute_expected_cis {
    tag "${label} (${resolution})"

    input:
        tuple val(label),
              val(resolution),
              path(cooler)

    output:
        tuple val(label),
              val(resolution),
              path("*.tsv"), emit: tsv

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        cooltools expected-cis \\
            "$cooler" \\
            -o '!{outprefix}_expected_cis.tsv'
        '''
}

process call_loops {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_medium'

    tag "${label} (${resolution})"

    input:
        tuple val(label),
              val(resolution),
              path(expected_cis),
              path(cooler)

    output:
        tuple val(label),
              val(resolution),
              path("*.bedpe.gz"), emit: bedpe

    shell:
        outprefix="${label}_${resolution}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        cooltools dots \\
            "$cooler" \\
            '!{expected_cis}' \\
            --max-loci-separation 6000000 \\
            --clustering-radius 50000 \\
            -p '!{task.cpus}' \\
            -o '!{label}'

        gzip -9c '!{label}' > '!{outprefix}.bedpe.gz'
        '''
}

process merge_loops {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    tag "${label}"

    input:
        tuple val(label),
              val(resolutions),
              path(loops)

        val lowest_resolution

    output:
        tuple val(label),
              path("*.bedpe.gz"), emit: bedpe

    shell:
        '''
        for f in *.gz; do
            # Convert bedpe to bedgraph and remove header
            zcat "$f" |
            cut -f 1-6 |
            tail -n +2 > "$f.input"
        done

        hicMergeLoops --inputFiles *.input \\
            --outFileName '!{label}' \\
            --lowestResolution '!{lowest_resolution}'

        gzip -9c '!{label}' > '!{label}.bedpe.gz'
        '''
}

process detect_differential_loops {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path loops
        val lowest_resolution

    output:
        path "*.tsv.gz", emit: tsv

    shell:
        '''
        identify_differential_loops.py \\
            *_merged.bedpe.gz \\
            '!{lowest_resolution}' |
            gzip -9 > differential_loops_by_condition.tsv.gz

        identify_differential_loops.py \\
            hg38_00*.bedpe.gz \\
            '!{lowest_resolution}' |
            gzip -9 > differential_loops_by_sample.tsv.gz
        '''

}

process generate_loop_report {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    label 'process_medium'

    input:
        path mcools
        path loops
        val lowest_resolution

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg

    shell:
        '''
        generate_loop_report.py \\
            --mcools *{WT,T1,C1}*_merged.mcool \\
            --loops *{WT,T1,C1}_merged.bedpe.gz \\
            --lowest-resolution '!{lowest_resolution}' \\
            -o loops \\
            --labels WT T1 C1
        '''
}
