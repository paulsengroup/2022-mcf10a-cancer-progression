#!/usr/bin/env nextflow
// Copyright (C) 2023 Saleh Oshaghi <mohao@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.translocations, checkIfExists: true)
           .map { tuple(it.getName().replaceAll(/_translocations.bed$/, ""), file(it)) }
           .set { translocations }

    Channel.fromPath(params.domains, checkIfExists: true)
           .map { tuple(it.getName().replaceAll(/_all_domains.bed.gz$/, ""), file(it)) }
           .set { domains }

    Channel.fromPath(params.cliques, checkIfExists: true)
           .map { tuple(it.getName().replaceAll(/_all_cliques.tsv.gz$/, ""), file(it)) }
           .set { cliques }

    Channel.fromPath(params.lads, checkIfExists: true)
           .map { tuple(it.getName().replaceAll(/_LMNB1_unionrep_peaks.bed.gz$/, ""), file(it)) }
           .set { lads }

    generate_blacklist(translocations,
                       file(params.blacklist, checkIfExists: true))

    domains\
        .join(cliques, failOnDuplicate: true, failOnMismatch: true)
        .join(lads, failOnDuplicate: true, failOnMismatch: true)
        .join(generate_blacklist.out.bed, failOnDuplicate: true, failOnMismatch: true)
        .set { preprocess_significant_interactions_input_ch }

    preprocess_significant_interactions(
        preprocess_significant_interactions_input_ch,
        file(params.chr_sizes, checkIfExists: true),
        params.bin_size)

    generate_seed_sequence(
        preprocess_significant_interactions.out.gtrack
            .map { it[1] }  // Discard label
            .collect(),
            params.number_of_models)

    generate_seed_sequence.out.txt
        .splitText()
        .map { tuple(it.trim().split("\t")) }
        .set { seeds }

    preprocess_significant_interactions.out.gtrack
        .combine(seeds)
        .set { chrom3d_input_ch }

    run_chrom3d(chrom3d_input_ch,
                params.N,
                params.L,
                params.r)
}


process generate_blacklist {
    label 'process_very_short'

    input:
        tuple val(label),
              path(translocations)
        path blacklist

    output:
        tuple val(label),
              path("blacklist.bed"),
        emit: bed

    shell:
        '''
        set -o pipefail

        cat <(zcat '!{blacklist}') \\
            <(zcat '!{translocations}') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 > blacklist.bed
        '''
}

process preprocess_significant_interactions{
    label 'process_very_short'

    input:
        tuple val(label),
              path(domains),
              path(cliques),
              path(lads),
              path(blacklist)
        path chr_sizes
        val bin_size

    output:
        tuple val(label),
              path("*.gtrack"),
        emit: gtrack

    shell:
        outprefix="${domains.simpleName}"
        '''
        generate_list_of_interacting_domains_from_cliques.py '!{domains}' '!{cliques}' > '!{outprefix}.bedpe'

        chrom3d_tad_to_gtrack.sh \\
            '!{outprefix}.bedpe' \\
            '!{bin_size}' \\
            '!{chr_sizes}' \\
            '!{lads}' \\
            '!{blacklist}' > '!{outprefix}_LADs_diploid.gtrack'
        '''
}

process generate_seed_sequence {
    label 'process_very_short'

    input:
        path files
        val num_seeds

    output:
        path "seeds.txt", emit: txt

    shell:
        '''
        generate_seed_sequence.py !{files} --number-of-seeds='!{num_seeds}' > seeds.txt
        '''
}


process run_chrom3d {
    publishDir "${params.outdir}/cmms", mode: 'copy'

    label 'error_retry'

    input:
        tuple val(label),
              path(input_gtrack),
              val(id),
              val(seed)
        val N
        val L
        val radius

    output:
        tuple val(label),
              path("*.cmm"),
        emit: cmm

        tuple val(label),
              val(seed),
        emit: seed

    shell:
        outname="${input_gtrack.simpleName}_${id}.cmm"
        '''
        Chrom3D -o '!{outname}' \\
                -r '!{radius}' \\
                -n '!{N}' \\
                -l '!{L}' \\
                -s '!{seed}' \\
                --nucleus '!{input_gtrack}'
        '''
}
