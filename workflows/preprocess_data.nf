#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    sort_and_filter_chrom_sizes(
        file(params.hg38_chrom_sizes_in, checkIfExists: true),
        params.hg38_chrom_sizes_out
    )

    sort_and_filter_fna(
        file(params.hg38_assembly_in, checkIfExists: true),
        params.hg38_assembly_filtered_out
    )

    process_microarray_data(
        file(params.hg38_chrom_sizes_in, checkIfExists: true),
        file(params.microarray_cnvs, checkIfExists: true),
        file(params.microarray_probe_dbs, checkIfExists: true),
        file(params.hg17_to_hg38_liftover_chain, checkIfExists: true)
    )

    generate_blacklist(
        file(params.hg38_assembly_gaps, checkIfExists: true),
        file(params.hg38_cytoband, checkIfExists: true),
        params.hg38_blacklist
    )

    Channel.of(
        tuple(file(params.hg38_assembly_in), params.hg38_assembly_out),
        tuple(file(params.hg38_gtf_in), params.hg38_gtf_out)
    ).set { files }

    decompress_data(
        files
    )

    untar_data(
        tuple(file(params.fish_tar), params.fish_out_dir)
    )
}

process sort_and_filter_chrom_sizes {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${chrom_sizes_out}" }

    label 'process_very_short'

    input:
        path chrom_sizes_in
        val chrom_sizes_out

    output:
        path "*.chrom.sizes", emit: chrom_sizes

    shell:
        out=file(chrom_sizes_out).getName()
        '''
        grep '^chr[[:digit:]XY]\\+[[:space:]]' '!{chrom_sizes_in}' |
           sort -V > '!{out}'
        '''
}

process sort_and_filter_fna {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${fna_out}" }

    input:
        path fna_in
        val fna_out

    output:
        path "*.fa", emit: fa

    shell:
        out=file(fna_out).getName()
        '''
        seqkit grep  \\
            --pattern '^chr[XY\\d]+$' \\
            --use-regexp \\
            '!{fna_in}' > fa.tmp

        seqkit sort \\
            --by-name \\
            --natural-order \\
            --two-pass \\
            fa.tmp > '!{out}'
        '''
}

process generate_blacklist {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${dest}" }

    label 'process_very_short'

    input:
        path assembly_gaps
        path cytoband
        val dest

    output:
        path "*.bed.gz", emit: bed

    shell:
        '''
        set -o pipefail

        cat <(zcat '!{assembly_gaps}' | cut -f 2-) \\
            <(zcat '!{cytoband}' | grep 'acen$') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 |
            gzip -9 > blacklist.bed.gz
        '''
}

process process_microarray_data {
    publishDir "${params.output_dir}/microarray", mode: 'copy'
    label 'process_short'

    input:
        path chrom_sizes
        path bed
        path probe_dbs
        path liftover_chain

    output:
        path "*.bed.gz", emit: bed

    shell:
        outname="${bed.simpleName}.bed.gz"
        '''
        set -o pipefail
        convert_microarray_cnvs_to_bed.py \\
            '!{bed}' \\
            --chrom-sizes '!{chrom_sizes}' \\
            --probe-ids !{probe_dbs} \\
            --liftover-chain '!{liftover_chain}' \\
            --fill-gaps |
            gzip -9c > '!{outname}'
        '''
}

process decompress_data {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${dest}" }

    label 'process_short'

    input:
        tuple path(src),
              val(dest)

    output:
        path "*", emit: file

    shell:
        out=file(dest).getName()
        '''
        zcat '!{src}' > '!{out}'
        '''
}


process untar_data {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${dest}/${it}" }

    label 'process_short'

    input:
        tuple path(src),
              val(dest)

    output:
        path "*", type: "dir", emit: file

    shell:
        '''
        tar -xf '!{src}'
        '''
}
