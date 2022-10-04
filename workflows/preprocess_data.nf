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

    run_bowtie2_index(rename_chromosomes.out.fa)
}

process generate_chrom_sizes {
    publishDir "${params.output_dir}", mode: 'copy'

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
         awk -F $'\t' 'BEGIN { OFS=FS } $2 == "assembled-molecule" { print "chr"$1,0,$9,$7 }' |
         grep -v 'chrMT' | sort -V > "!{out_bed}"

         # Convert chromosome sizes from bed to chrom.sizes format
         cut -f 1,3 "!{out_bed}" | sort -V > "!{out}"
        '''
}

process run_bowtie2_index {
    publishDir "${params.output_dir}/bowtie2_idx/", mode: 'copy'

    label 'process_high'

    input:
        path fa

    output:
        path "*.bt2"

    shell:
        '''
        fa='!{fa}'
        outprefix="${fa%.*}"

        bowtie2-build --threads !{task.cpus} \
                      "$fa"                  \
                      "$outprefix"
        '''
}

process rename_chromosomes {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path fa
        path chrom_sizes_bed

    output:
        path "${fa.baseName}", emit: fa

    shell:
        out="${fa.baseName}"
        '''
        gzip -dc '!{fa}' | '!{params.script_dir}/rename_chromosomes_fa.py' '!{chrom_sizes_bed}' > '!{out}'
        '''
}
