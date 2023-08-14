#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


def collect_files(prefix, sample_id, suffix, type = "file") {
    // Example: MCF10A_WT_REP3_T1 -> MCF10A_WT
    def files = file("${prefix}/${sample_id}*${suffix}",
                     type: type,
                     checkIfExists: true)

    def condition_id = sample_id.replaceAll(/(.*)_REP\d(.*)/, '$1')


    return tuple(sample_id, condition_id, files)
}


workflow {
    Channel.fromPath(params.nfcore_samplesheet)
        .splitCsv(sep: ",", header: true)
        .map { it.sample.replaceAll(/_T1$/, '') }
        .unique()
        .set { sample_ids }

    input_dir = file(params.input_dir, type: 'dir', checkIfExists: true)

    sample_ids
        .map { collect_files("${input_dir}/star_salmon/bigwig/", it, '.bigWig') }
        .set { input_bigwigs }

    reorder_chromosomes_bigwig(
        input_bigwigs,
        file(params.chrom_sizes, checkIfExists: true)
    )

    merge_bigwigs(
        reorder_chromosomes_bigwig.out.bwig
            .groupTuple(by: 1)
                        // condition, bigwigs
            .map{ tuple(it[1], it[2].flatten()) },
        file(params.chrom_sizes, checkIfExists: true)
    )

}

process reorder_chromosomes_bigwig {
    publishDir "${params.output_dir}/star_salmon/bigwig/sorted", mode: 'copy'
    tag "$sample"

    input:
        tuple val(sample),
              val(condition),
              path(bwigs)
        path chrom_sizes

    output:
        tuple val(sample),
              val(condition),
              path("*.sorted.bw"),
        emit: bwig

    shell:
        '''
        for bwin in !{bwigs}; do
            bwout="${bwin%.bigWig}"
            bwout="${bwin%.bigwig}"
            bwout="${bwout%.bw}.sorted.bw"

            reorder_chromosomes_bigwig.py \\
                "$bwin" \\
                '!{chrom_sizes}' \\
                -o "$bwout"
        done
        '''
}

process merge_bigwigs {
    publishDir "${params.output_dir}/star_salmon/bigwig/merged", mode: 'copy'
    tag "$condition"

    input:
        tuple val(condition),
              path(bwigs)
        path chrom_sizes

    output:
        tuple val(condition),
              path("*.bw"),
        emit: bwig

    shell:
        outname="${condition}.bw"
        '''
        wiggletools sum !{bwigs} > tmp.wig
        wigToBigWig tmp.wig '!{chrom_sizes}' '!{outname}'
        '''
}
