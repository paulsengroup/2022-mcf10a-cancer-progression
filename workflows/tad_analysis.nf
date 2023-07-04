#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.mcools)
           .map {
                if ("${it}".endsWith("_merged.mcool")) {
                    return tuple("cond", it)
                }
                return tuple("sample", it)
           }
           .set{ mcools }

    normalization_methods = params.norm_methods.split(",")
    resolutions = "${params.resolutions}".split(",")

    apply_normalization_to_coolers(
        mcools.combine(normalization_methods)
              .combine(Channel.of(resolutions))
              .map { tuple(it[0], it[2], it[3], it[1]) },
              // sample_type, norm, resolution, cooler
        file(params.blacklist, checkIfExists: true))

    hicexplorer_find_tads(apply_normalization_to_coolers.out.cool)
    hicexplorer_find_tads.out.tads
                         .groupTuple(by: [1, 2])
                         .set { tads }

    generate_chrom_sizes(
        mcools.first(),
        resolutions.first()
    )

    insulation_to_bigwig(
        hicexplorer_find_tads.out.tads,
        generate_chrom_sizes.out.chrom_sizes
    )

    aggregate_insulation_scores(
        tads,
        params.repl_labels,
        params.cond_labels
    )

    generate_tad_report(
        tads,
        params.repl_pretty_labels,
        params.cond_pretty_labels
    )

    generate_insulation_report(
        aggregate_insulation_scores.out.bedgraph,
        tads.map { it[5] }.collect(), // domains
        file(params.blacklist),
        params.repl_pretty_labels,
        params.cond_pretty_labels
    )
}

process generate_chrom_sizes {

    input:
        tuple val(sample_type), path(cooler)
        val resolution

    output:
        path "*.chrom.sizes", emit: chrom_sizes

    shell:
        outprefix="${cooler.baseName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi
        cooler dump "$cooler" -t chroms > '!{outprefix}.chrom.sizes'
        '''
}

process apply_normalization_to_coolers {
    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(cooler)
        path blacklist

    output:
        tuple val(type),
              val(normalization),
              val(resolution),
              path("*.cool"), emit: cool

    shell:
        out="${cooler.baseName}_${normalization}_${resolution}.cool"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        cooler_apply_normalization.py              \\
                    --norm-name='!{normalization}' \\
                    --output-name='!{out}'         \\
                    --nproc='!{task.cpus}'         \\
                    "$cooler"
        '''
}

process hicexplorer_find_tads {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: {
                                             it = it.toString().replaceAll("_${normalization}_${resolution}_", "_")
                                            "${normalization}/${resolution}/${it}"
                                       }

    label 'process_medium'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(cooler)

    output:
        tuple val(type),
              val(normalization),
              val(resolution),
              path("*boundaries.bed.gz"),
              path("*domains.bed.gz"),
              path("*score.bedgraph.gz"),
              path("*tad_score.bm.gz"), emit: tads

    shell:
        outprefix="${cooler.baseName}"
        '''
        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        hicFindTADs -p '!{task.cpus}'          \\
                    --matrix '!{cooler}'       \\
                    --outPrefix '!{outprefix}' \\
                    --correctForMultipleTesting fdr

        # Compress text files
        printf '%s\\n' *.bed *.bedgraph *.bm |
            xargs -L 1 -P '!{task.cpus}' sh -c 'gzip -9 "$1"' sh
        '''
}

process insulation_to_bigwig {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: {
                                             it = it.toString().replaceAll("_${normalization}_${resolution}_", "_")
                                            "${normalization}/${resolution}/${it}"
                                       }

    label 'process_very_short'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(boundaries),
              path(domains),
              path(insulation_score),
              path(tad_score)
        path chrom_sizes

    output:
        path "*.bw", emit: bw

    shell:
        outprefix="${insulation_score.simpleName}"
        '''
        zcat '!{insulation_score}' |
        sort -k1,1 -k2,2n > scores.bedgraph

        bedGraphToBigWig scores.bedgraph '!{chrom_sizes}' '!{outprefix}.bw'
        '''
}

process aggregate_insulation_scores {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    label 'process_short'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(boundaries),
              path(domains),
              path(insulation_score),
              path(tad_score)

        val labels_replicates
        val labels_conditions

    output:
        tuple val(type),
              val(normalization),
              val(resolution),
              path("*.bedgraph.gz"), emit: bedgraph

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        aggregate_insulation_scores.py \\
            *WT*merged!{suffix}_domains.bed.gz \\
            --labels='!{labels_replicates}' |
             gzip -9 > 'insulation_scores_replicates.bedgraph.gz'

        aggregate_insulation_scores.py \\
            *WT*merged!{suffix}_domains.bed.gz \\
            --labels='!{labels_conditions}' |
             gzip -9 > 'insulation_scores_conditions.bedgraph.gz'
        '''
}


process generate_tad_report {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    label 'process_short'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(boundaries),
              path(domains),
              path(insulation_scores),
              path(tad_score)

        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        generate_tad_report.py \
            hg38_???_*!{suffix}_domains.bed.gz \
            --output-prefix=report_replicates \
            --labels='!{labels_replicates}'

        generate_tad_report.py \
            *{WT,T1,C1}_*merged!{suffix}_domains.bed.gz \
            --output-prefix=report_conditions \
            --labels='!{labels_conditions}'

        mkdir plots
        mv *.svg *.png plots/
        '''
}

process generate_insulation_report {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    label 'process_short'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(insulation)

        path tads
        path blacklist

        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        if [[ '!{type}' == cond ]]; then
            labels='!{labels_conditions}'
        else
            labels='!{labels_replicates}'
        fi

        generate_insulation_report.py \\
            '!{insulation}' \\
            *WT_*merged!{suffix}_domains.bed.gz \
            --output-prefix=report_replicates_insulation \\
            --blacklist '!{blacklist}' \\
            --labels="$labels"

        mkdir plots
        mv *.svg *.png plots/
        '''
}
