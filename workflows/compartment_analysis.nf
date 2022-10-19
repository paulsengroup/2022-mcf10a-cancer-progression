#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    call_compartments(Channel.fromPath(params.mcools),
                      file(params.ref_genome),
                      params.resolution)
}

process call_compartments {
    publishDir "${params.output_dir}/compartments", mode: 'copy'

    input:
        path cooler
        path ref_genome
        val resolution

    output:
        path "*.vecs.tsv", emit: eigvect_txt
        path "*.lam.txt", emit: eigval_txt
        path "*.bw", emit: bw

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [[ '!{ref_genome}' == *.gz ]]; then
            zcat '!{ref_genome}' > ref.fna
        else
            ln -s '!{ref_genome}' ref.fna
        fi

        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        # Compute GC content. Used to orient eigvects
        '!{params.script_dir}/compute_binned_gc_content.py' \
            <(cooler dump -t bins "$cooler")                \
            '!{ref_genome}' > gc.bed

        cooltools eigs-cis --phasing-track gc.bed \
                           --bigwig               \
                           -o '!{outprefix}'      \
                           "$cooler"
        '''
}

// TODO: generate saddleplots: https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-saddle
