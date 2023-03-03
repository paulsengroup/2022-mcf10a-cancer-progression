#!/usr/bin/env nextflow
// Copyright (C) 2023 Saleh Oshaghi <mohao@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    preprocess_significant_interactions(Channel.fromPath(params.domains),
    Channel.fromPath(params.cliques),
    file(params.lads), file(params.cytoband), file(params.gaps), file(params.chr_sizes), params.bin_size, file(params.translocation), params.outdir)

    run_chrom3d(preprocess_significant_interactions.out.gtrack, params.N, params.L, params.r)
}

process preprocess_significant_interactions{
    publishDir "${params.outdir}", mode: 'copy'
    input:
        //tuple val(sampleId), path(read)
        path domain
        path clique
        path lads
        path cytoband
        path gaps
        path chr_sizes
        val bin_size
        path translocation
        path outdir
    output:
        path "*.gtrack", emit: gtrack
    shell:
    '''
    generate_cliques_bedpe.py !{domain} !{clique} > T1_output_significant.txt
    chrom3d_tad_to_gtrack.sh T1_output_significant.txt '!{bin_size}' '!{chr_sizes}' '!{cytoband}' '!{gaps}' '!{lads}' '!{translocation}'
    '''
}


process run_chrom3d{
    publishDir "${params.outdir}", mode: 'copy'
    input:
        path input_gtrack
        val N
        val L
        val radius
    output:
        path "*.cmm"
    shell:
    '''
    Chrom3D -o !{input_gtrack}.cmm -r !{radius} -n !{N} -l !{L} --nucleus !{input_gtrack}
    '''
}
