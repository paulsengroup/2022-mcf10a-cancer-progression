#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.fromList(params.mcools_by_condition)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { coolers_by_condition }

    Channel.fromList(params.mcools_by_sample)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { coolers_by_sample }

    Channel.fromList(params.hint_cnv_profiles)
        .map { tuple(it[0], file(it[1], checkIfExists: true)) }
        .set { hint_cnv_tars }

    run_hicrep(
        coolers_by_sample
            .map { it[1] }
            .collect(),
        params.plot_pretty_labels,
        params.hicrep_bin_size,
        params.hicrep_h
    )

    extract_chrom_sizes_from_mcool(
        coolers_by_sample.first(),
        params.hicrep_bin_size
    )

    extract_chrom_sizes_from_mcool.out.chrom_sizes
        .set { chrom_sizes }

    plot_scc(
        chrom_sizes,
        run_hicrep.out.tsv
    )

    run_maxhic(
        coolers_by_condition,
        params.maxhic_bin_size
    )

    merge_hint_cnv_profiles(
        hint_cnv_tars,
        chrom_sizes
    )

    coolers_by_condition
        .join(merge_hint_cnv_profiles.out.bedgraph)
        .set { balancing_tasks }


    balance_with_loic(
        balancing_tasks,
        params.cnv_balancing_resolution,
        params.cnv_balancing_chromosomes
    )

    balance_with_caic(
        balancing_tasks,
        params.cnv_balancing_resolution,
        params.cnv_balancing_chromosomes
    )
}

process extract_chrom_sizes_from_mcool {
    label 'process_very_short'

    input:
        tuple val(sample),
              path(mcool)

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

        run_hicrep.py \\
            "${coolers[@]}" \\
            -p '!{task.cpus}' \\
            --h '!{h}' \\
            --labels='!{labels}' |
            grep -v 'chr[MY]' > hicrep.scc.tsv
        '''
}

process plot_scc {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    label 'process_very_short'

    input:
        path chrom_sizes
        path hicrep_output

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg

    shell:
        '''
        plot_scc.py \\
            '!{hicrep_output}' \\
            --chrom-sizes '!{chrom_sizes}' \\
            -o hicrep
        '''
}

process run_maxhic {
    publishDir "${params.output_dir}/maxhic", mode: 'copy'

    label 'process_high'

    input:
        tuple val(sample),
              path(mcool)

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

process merge_hint_cnv_profiles {
    label 'process_short'

    input:
        tuple val(sample),
              path(tar)

        path chrom_sizes

    output:
        tuple val(sample),
              path(outname),
        emit: bedgraph

    shell:
        outname="${sample}.hint.cnv.bedgraph.gz"
        '''
        set -o pipefail

        mkdir hint-out
        tar -xf '!{tar}' -C hint-out --strip-components=1

        merge_hint_cnv_profiles.py \\
            hint-out/segmentation/b*/ \\
            '!{chrom_sizes}' |
            gzip -9 > '!{outname}'
        '''
}

process balance_with_loic {
    publishDir "${params.output_dir}/balanced_matrices/loic/", mode: 'copy'

    label 'process_long'

    input:
        tuple val(sample),
              path(mcool),
              path(cnvs)

        val resolution
        val chroms

    output:
        tuple val(sample),
              path(outname),
        emit: cool

    shell:
        uri=mcool
        if(resolution != 0) {
            uri="${mcool}::/resolutions/${resolution}"
        }
        outname="${sample}.${resolution}.loic.cool"
        '''
        cooler_balance_cnv.py \\
            '!{uri}' \\
            '!{cnvs}' \\
            '!{outname}' \\
            --method=LOIC \\
            --chromosomes='!{chroms}'
        '''
}

process balance_with_caic {
    publishDir "${params.output_dir}/balanced_matrices/caic/", mode: 'copy'

    label 'process_long'

    input:
        tuple val(sample),
              path(mcool),
              path(cnvs)

        val resolution
        val chroms

    output:
        tuple val(sample),
              path(outname),
        emit: cool

    shell:
        uri=mcool
        if(resolution != 0) {
            uri="${mcool}::/resolutions/${resolution}"
        }
        outname="${sample}.${resolution}.caic.cool"
        '''
        cooler_balance_cnv.py \\
            '!{uri}' \\
            '!{cnvs}' \\
            '!{outname}' \\
            --method=CAIC \\
            --chromosomes='!{chroms}'
        '''
}
