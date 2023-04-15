#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    grch38_bname = "${params.grch38_assembly_name_short}"

    bnames = channel.of(grch38_bname)

    assembly_reports = channel.of(file(params.grch38_assembly_report))
    genome_assemblies = channel.of(file(params.grch38_genome_assembly))

    generate_chrom_sizes(bnames,
                         assembly_reports)

    rename_chromosomes(genome_assemblies,
                       generate_chrom_sizes.out.bed)

    process_microarray_data(generate_chrom_sizes.out.chrom_sizes,
                            file(params.microarray_cnvs),
                            file(params.microarray_probe_dbs),
                            file(params.hg17_to_hg38_liftover_chain))

    decompress_data(
        Channel.of(
            file(params.hg38_assembly_in, checkIfExists: true),
            file(params.hg38_gtf_in, checkIfExists: true)
        ),
        Channel.of(
            params.hg38_assembly_out,
            params.hg38_gtf_out
        )
    )

    split_fastq_pair(
        Channel.fromFilePairs(params.hic_fastq_pattern,
                              checkIfExists: true,
                              flat: true),
        params.fastq_chunk_size)
}

process generate_chrom_sizes {
    publishDir "${params.output_dir}/chrom_sizes", mode: 'copy'

    label 'process_short'

    input:
        val assembly_name
        path assembly_report

    output:
        path "${assembly_name}.bed", emit: bed
        path "${assembly_name}.chrom.sizes", emit: chrom_sizes

    shell:
        out_bed = "${assembly_name}.bed"
        out = "${assembly_name}.chrom.sizes"
        '''
        set -e
        set -u
        set -o pipefail

         # Extract chromosome sizes from assembly report
         gzip -dc "!{assembly_report}" |
         awk -F $'\\t' 'BEGIN { OFS=FS } $2 == "assembled-molecule" { print "chr"$1,0,$9,$7 }' |
         grep -v 'chrMT' | sort -V > "!{out_bed}"

         # Convert chromosome sizes from bed to chrom.sizes format
         cut -f 1,3 "!{out_bed}" | sort -V > "!{out}"
        '''
}


process rename_chromosomes {
    publishDir "${params.output_dir}/assemblies", mode: 'copy'

    label 'process_short'

    input:
        path fa
        path chrom_sizes_bed

    output:
        path "${fa.baseName}", emit: fa

    shell:
        out="${fa.baseName}"
        '''
        gzip -dc '!{fa}' | rename_chromosomes_fa.py '!{chrom_sizes_bed}' > '!{out}'
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
    publishDir "${params.output_dir}", mode: 'copy',
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

process split_fastq_pair {
    publishDir "${params.output_dir}/fastq", mode: 'copy'
    label 'process_long'

    cpus 3

    input:
        tuple val(key),
              path(m1),
              path(m2)
        val chunk_size

    output:
        path "*.fastq.zst"

    shell:
    '''
    seqkit split2 -s '!{chunk_size}' \\
        -1 '!{m1}' \\
        -2 '!{m2}' \\
        -e ".zst"
    '''
}
