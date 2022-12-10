#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    coolers_by_condition = Channel.fromPath(params.mcools_by_condition)
    coolers_by_sample = Channel.fromPath(params.mcools_by_sample)

    run_hicrep(coolers_by_sample.collect(),
               params.plot_pretty_labels,
               params.hicrep_bin_size,
               params.hicrep_h)

    extract_chrom_sizes_from_mcool(coolers_by_sample.first(), params.hicrep_bin_size)

    plot_scc(extract_chrom_sizes_from_mcool.out.chrom_sizes,
             run_hicrep.out.tsv)

    run_maxhic(coolers_by_condition,
               params.maxhic_bin_size)
}

process extract_chrom_sizes_from_mcool {
    label 'very_short'

    input:
        path mcool
        val resolution

    output:
        path "chrom.sizes", emit: chrom_sizes

    shell:
        '''
        cooler dump -t chroms '!{mcool}::/resolutions/!{resolution}' > chrom.sizes
        '''
}

process run_hicrep {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_high'

    input:
        path coolers
        val labels
        val resolution
        val h

    output:
        path "*.tsv", emit: tsv

    shell:
        '''
        set -o pipefail
        shopt -s nullglob

        coolers=()
        if [ '!{resolution}' -ne 0 ]; then
            for cooler in *.mcool *.cool; do
                coolers+=("$cooler::/resolutions/!{resolution}")
            done
        else
            coolers=(!{coolers})
        fi

        '!{params.script_dir}/run_hicrep.py' \
            "${coolers[@]}" \
            -p '!{task.cpus}' \
            --h '!{h}' \
            --labels='!{labels}' |
            grep -v 'chr[MY]' > hicrep.scc.tsv
        '''
}

process plot_scc {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    label 'very_short'

    input:
        path chrom_sizes
        path hicrep_output

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg

    shell:
        '''
        '!{params.script_dir}/plot_scc.py' \
            '!{hicrep_output}' \
            --chrom-sizes '!{chrom_sizes}' \
            -o hicrep
        '''
}

process run_maxhic {
    publishDir "${params.output_dir}/maxhic", mode: 'copy'

    label 'process_high'

    input:
        path mcool
        val resolution

    output:
        path "*cis_interactions.tsv.zst", emit: cis_interactions
        path "*trans_interactions.tsv.zst", emit: trans_interactions
        path "*model_params.tar.zst", emit: model_params

    shell:
        outprefix="${mcool.baseName}"
        '''
        set -o pipefail

        trap 'rm -rf input/ output/' EXIT
        mkdir input output/

        cooler dump -t bins '!{mcool}::/resolutions/!{resolution}' |
            awk -F '\\t' 'BEGIN{ OFS = FS } {print $1,$2,$3,NR}' > input/bins.bed

        cooler dump -t pixels --one-based-ids '!{mcool}::/resolutions/!{resolution}' > input/pixels.matrix

        maxhic -t !{task.cpus} input/ output/

        zstd -T!{task.cpus} -13 output/cis_interactions.txt -o '!{outprefix}_cis_interactions.tsv.zst'
        zstd -T!{task.cpus} -13 output/trans_interactions.txt -o '!{outprefix}_trans_interactions.tsv.zst'

        mv output '!{outprefix}'

        tar -cf - '!{outprefix}/ModelParameters/' | zstd -T!{task.cpus} -13 -o '!{outprefix}_model_params.tar.zst'
        '''
}
