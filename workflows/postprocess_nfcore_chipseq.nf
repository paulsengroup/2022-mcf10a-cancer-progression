#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


def collect_files(prefix, sample_id, suffix, type = "file") {
    def files = file("${prefix}/${sample_id}${suffix}",
                     type: type,
                     checkIfExists: true)
    return tuple(sample_id, files)
}


workflow {
    Channel.fromPath(params.nfcore_samplesheet)
        .splitCsv(sep: ",", header: true)
        .filter { it.control != "" }
        .map {
            tuple(it.sample,
                  it.control)
             }
        .unique()
        .map {
            tuple(
                collect_files("${params.input_dir}/*/mergedLibrary/", it[0], '.mLb.clN.sorted.bam'),
                collect_files("${params.input_dir}/*/mergedLibrary/", it[1], '.mLb.clN.sorted.bam')
            )
            }
        .map { tuple(it[0][0], it[0][1], it[1][0], it[1][1]) }
        .set { input_bams }

    macs2_fold_change_over_control(
        input_bams
    )

    bedgraph_to_bigwig(
        macs2_fold_change_over_control.out.bedgraph
            .map { it[2] },
        file(params.chrom_sizes, checkIfExists: true)
    )
}

process macs2_fold_change_over_control {
    tag "${control}_${treatment}"

    input:
        tuple val(treatment),
              path(treatment_bam),
              val(control),
              path(control_bam)

    output:
        tuple val(treatment),
              val(control),
              path("*_fc_over_control.bdg"),
        emit: bedgraph

    shell:
        '''
        # See
        # https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/ec4295c8ac68be25b25357038d82ec942ac0bf8d/src/encode_task_macs2_signal_track_chip.py#L90C16-L102

        macs2 callpeak \\
            -g hs \\
            -t '!{treatment_bam}' \\
            -c '!{control_bam}' \\
            -n out \\
            -p 0.01 \\
            --nomodel \\
            --shift 0 \\
            --extsize 155 \\
            --keep-dup all \\
            --bdg \\
            --SPMR

        # See
        # https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/ec4295c8ac68be25b25357038d82ec942ac0bf8d/src/encode_task_macs2_signal_track_chip.py#L147C1-L154C6
        macs2 bdgcmp \\
            -t out_treat_pileup.bdg \\
            -c out_control_lambda.bdg \\
            --o-prefix out \\
            -m ppois

        mkdir tmp
        sort -T tmp/ \\
            -k1,1 -k2,2n \\
            out_ppois.bdg > '!{treatment}_fc_over_control.bdg'
        '''
}


process bedgraph_to_bigwig {
    publishDir "${params.output_dir}/chromap/mergedLibrary/bigwig/fold_change", mode: 'copy'

    tag "${bedgraph.simpleName}"

    input:
        path bedgraph
        path chrom_sizes

    output:
        path "*.bw", emit: bigwig

    shell:
        outname="${bedgraph.simpleName}.bw"
        '''
        grep -P 'chr[XY\\d+]\\s' '!{bedgraph}' > bedgraph.tmp

        bedGraphToBigWig \\
            bedgraph.tmp \\
            '!{chrom_sizes}' \\
            '!{outname}'
        '''
}
