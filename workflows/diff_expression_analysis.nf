#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.count_matrices, checkIfExists: true)
        .map {
            def bname = file(it).getName()
            def dirname = file(it).getParent()
            def count_type = bname =~ /(gene|transcript)_counts/
            def pipeline_type = file(dirname).getName()

            return tuple(count_type[0][1], pipeline_type, file(it))
        }.set { count_tables }

    count_tables
        .combine(Channel.of(params.lfc_cutoffs).flatten())
        .map {
            tuple(it[0],
                  it[1],
                  it[2],
                  file(params.design_table),
                  it[3])
        }.set { run_deseq2_input_ch }

    run_deseq2(run_deseq2_input_ch)
}


process run_deseq2 {
    publishDir "${params.output_dir}/", mode: 'copy'

    input:
        tuple val(count_type),
              val(pipeline_type),
              path(count_table),
              path(design_table),
              val(lfc_cutoff)

    output:
        tuple val(count_type),
              val(pipeline_type),
              val(lfc_cutoff),
              path("${pipeline_type}_${count_type}/lfc_${lfc_cutoff}/*.tsv.gz"), emit: tsv
        tuple val(count_type),
              val(pipeline_type),
              val(lfc_cutoff),
              path("${pipeline_type}_${count_type}/lfc_${lfc_cutoff}/*.rds"), emit: res
        tuple val(count_type),
              val(pipeline_type),
              val(lfc_cutoff),
              path("${pipeline_type}_${count_type}/lfc_${lfc_cutoff}/dds.rds"), emit: dds
        tuple val(count_type),
              val(pipeline_type),
              val(lfc_cutoff),
              path("${pipeline_type}_${count_type}/lfc_${lfc_cutoff}/r_sessioninfo.txt"), emit: sessioninfo

    shell:
        outdir = "${pipeline_type}_${count_type}/lfc_${lfc_cutoff}"
        '''
        run_deseq2.py          \\
            '!{count_table}'   \\
            '!{design_table}'  \\
            '!{outdir}'        \\
            --lfc-thresh='!{lfc_cutoff}'
        '''
}
