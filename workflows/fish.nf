#!/usr/bin/env nextflow
// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


workflow {

    Channel.fromList(params.images)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { images }

    Channel.fromList(params.num_green_blobs)
        .map { tuple(it[0], it[1]) }
        .set { num_green_blobs }

    Channel.fromList(params.num_red_blobs)
        .map { tuple(it[0], it[1]) }
        .set { num_red_blobs }


    segment_nuclei(
        images,
        params.contour_padding
    )

    segment_nuclei.out.h5
        .join(num_green_blobs)
        .join(num_red_blobs)
        .set { compute_blob_coordinates_tasks }

    compute_blob_coordinates(
        compute_blob_coordinates_tasks
    )

    plot_nuclei(
        compute_blob_coordinates.out.h5,
        params.overlap_cfx_cutoff
    )

    fetch_blob_coordinates(
        compute_blob_coordinates.out.h5,
        params.overlap_cfx_cutoff
    )
}


process segment_nuclei {
    input:
        tuple val(sample),
              path(pngs)

        val contour_padding

    output:
        tuple val(sample),
              path(outname),
        emit: h5

    shell:
    outname="${sample}_segmented_nuclei.h5"
    '''
    segment_nuclei.py \\
        *.png \\
        -o '!{outname}' \\
        --contour-padding='!{contour_padding}'
    '''
}

process compute_blob_coordinates {
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple val(sample),
              path(h5),
              val(num_green_blobs),
              val(num_red_blobs)

    output:
        tuple val(sample),
              path(outname),
        emit: h5

    shell:
    outname="${sample}_segmented_nuclei_with_blobs.h5"
    '''
    compute_blob_coordinates.py \\
        '!{h5}' \\
        '!{outname}' \\
        --num-green-blobs='!{num_green_blobs}' \\
        --num-red-blobs='!{num_red_blobs}'
    '''
}

process plot_nuclei {
    publishDir "${params.output_dir}/plots/", mode: 'copy'

    input:
        tuple val(sample),
              path(h5)

        val overlap_cfx

    output:
        tuple val(sample),
              path("*.pdf"),
        emit: pdf

    shell:
    outname1="${sample}.pdf"
    outname2="${sample}.blobs.filtered.pdf"
    '''
    plot_nuclei.py '!{h5}' -o '!{outname1}'

    plot_nuclei.py '!{h5}' \\
        -o '!{outname2}' \\
        --overlap-cfx-cutoff='!{overlap_cfx}' \\
        --plot-blobs
    '''
}

process fetch_blob_coordinates {
    publishDir "${params.output_dir}/blobs/", mode: 'copy'

    input:
        tuple val(sample),
              path(h5)

        val overlap_cfx

    output:
        tuple val(sample),
              path(outname),
        emit: tsv

    shell:
    outname="${sample}.blobs.tsv.gz"
    '''
    set -o pipefail

    fetch_blob_coordinates.py '!{h5}' --overlap-cfx-cutoff='!{overlap_cfx}' |
        gzip -9 > '!{outname}'
    '''
}
