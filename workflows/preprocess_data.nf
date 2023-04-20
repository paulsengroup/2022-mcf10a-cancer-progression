#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    filter_chrom_sizes(file(params.hg38_chrom_sizes_in, checkIfExists: true),
                       file(params.hg38_chrom_sizes_out, checkIfExists: true))

    run_bowtie2_index(rename_chromosomes.out.fa)
    archive_bowtie2_index(run_bowtie2_index.out.idx)

    process_microarray_data(file(params.hg38_chrom_sizes_in, checkIfExists: true),
                            file(params.microarray_cnvs, checkIfExists: true),
                            file(params.microarray_probe_dbs, checkIfExists: true),
                            file(params.hg17_to_hg38_liftover_chain, checkIfExists: true))

    decompress_data(
        Channel.fromPath([params.hg38_assembly_in, params.hg38_gtf_in], checkIfExists: true),
        Channel.fromPath([params.hg38_assembly_out, params.hg38_gtf_out], checkIfExists: false)
    )
}

process filter_chrom_sizes {
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
        grep '^chr[[:digit:]XY]\\+[[:space:]]' '!{chrom_sizes_in}' > '!{out}'
        '''
}

process run_bowtie2_index {
    label 'process_high'

    input:
        path fa

    output:
        path "*.bt2", emit: idx

    shell:
        '''
        fa='!{fa}'
        outprefix="${fa%.*}"

        bowtie2-build --threads !{task.cpus} \
                      "$fa"                  \
                      "$outprefix"
        '''
}

process archive_bowtie2_index {
    publishDir "${params.output_dir}/bowtie2_idx/", mode: 'copy'

    label 'process_medium'
    label 'process_short'

    input:
        path idx_files

    output:
        path "*.tar.zst", emit: tar

    shell:
        outname = "${idx_files[0].simpleName}.tar.zst"
        '''
        tar -chf - *.bt2 | zstd -T!{task.cpus} --adapt -o '!{outname}'
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
        convert_microarray_cnvs_to_bed.py \
            '!{bed}' \
            --chrom-sizes '!{chrom_sizes}' \
            --probe-ids !{probe_dbs} \
            --liftover-chain '!{liftover_chain}' \
            --fill-gaps |
            gzip -9c > '!{outname}'
        '''
}

process decompress_data {
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { "${dest}" }

    label 'process_short'

    input:
        path src
        val dest

    output:
        path "*", emit: file

    shell:
        out=file(dest).getName()
        '''
        zcat '!{src}' > '!{out}'
        '''
}
