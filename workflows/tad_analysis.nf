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
    hicexplorer_find_tads.out.tads  // normalization, resolution, tads
                         .map { tuple(it[1], it[2], it[4]) }
                         .groupTuple(by: [0, 1])
                         .set { tads }

    hicexplorer_find_tads.out.tads  // normalization, resolution, scores
                         .map { tuple(it[1], it[2], it[5]) }
                         .groupTuple(by: [0, 1])
                         .set { insulation }

    generate_chrom_sizes(
        mcools.first(),
        resolutions.first()
    )

    generate_blacklist(
        Channel.fromPath(params.mcools).collect(),
        Channel.of(normalization_methods)
            .combine(Channel.of(resolutions))
    )

    insulation_to_bigwig(
        hicexplorer_find_tads.out.tads,
        generate_chrom_sizes.out.chrom_sizes
    )

    aggregate_insulation_scores(
        insulation.join(generate_blacklist.out.bed, by: [0, 1]),
        params.repl_labels,
        params.cond_labels
    )

    generate_tad_report(
        tads,
        params.repl_pretty_labels,
        params.cond_pretty_labels
    )

    generate_tad_interaction_scatter(
        apply_normalization_to_coolers.out.cool
            .map { tuple(it[1], it[2], it[3]) }
            .groupTuple(by: [0, 1])
            .join(tads, by: [0, 1]),
        params.repl_labels,
        params.cond_labels
    )

    generate_insulation_report(
        aggregate_insulation_scores.out.bedgraph
            .join(generate_blacklist.out.bed, by: [0, 1]),
        tads.map { it[2] }.collect() // domains
    )

    generate_tad_overlap_report(
        tads.join(generate_blacklist.out.bed, by: [0, 1]),
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

process generate_blacklist {

    tag "${normalization}_${resolution}"

    input:
        path coolers
        tuple val(normalization),
              val(resolution)

    output:
        tuple val(normalization),
              val(resolution),
              path("*.bed.gz"),
        emit: bed

    shell:
        outname="${normalization}_${resolution}_blacklist.bed.gz"
        '''
        set -o pipefail

        uris=()
        for clr in !{coolers}; do
            if [ '!{resolution}' -ne 0 ]; then
                clr="$clr::/resolutions/!{resolution}"
            fi
            uris+=("$clr")
        done

        generate_blacklist_from_coolers.py ${uris[*]} |
            gzip -9 > '!{outname}'

        '''
}

process apply_normalization_to_coolers {
    tag "${cooler}_${normalization}_${resolution}"

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
    tag "${cooler}_${normalization}_${resolution}"

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
              path("*tad_score.bm.gz"),
              path("*tad_interactions.bedgraph.gz"),
        emit: tads

    shell:
        outprefix="${cooler.baseName}"
        '''
        NUMEXPR_MAX_THREADS='!{task.cpus}'
        export NUMEXPR_MAX_THREADS

        hicFindTADs -p '!{task.cpus}'          \\
                    --matrix '!{cooler}'       \\
                    --outPrefix '!{outprefix}' \\
                    --correctForMultipleTesting fdr

        annotate_tads_with_interactions.py '!{cooler}' \\
            *domains.bed.gz |
            gzip -9 > '!{outprefix}_tad_interactions.bedgraph.gz'

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
    tag "${normalization}_${resolution}"

    label 'process_very_short'

    input:
        tuple val(type),
              val(normalization),
              val(resolution),
              path(boundaries),
              path(domains),
              path(insulation_score),
              path(tad_score),
              path(interactions)
        path chrom_sizes

    output:
        path "*.bw", emit: bw

    shell:
        outprefix="${insulation_score.simpleName}"
        '''
        set -o pipefail

        zcat '!{insulation_score}' |
        sort -k1,1 -k2,2n > scores.bedgraph

        bedGraphToBigWig scores.bedgraph '!{chrom_sizes}' '!{outprefix}.bw'
        '''
}

process aggregate_insulation_scores {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }
    tag "${normalization}_${resolution}"

    label 'process_short'

    input:
        tuple val(normalization),
              val(resolution),
              path(insulation_scores),
              path(blacklist)

        val labels_replicates
        val labels_conditions

    output:
        tuple val(normalization),
              val(resolution),
              path("*.bedgraph.gz"), emit: bedgraph

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        set -o pipefail

        aggregate_insulation_scores.py \\
            hg38_???_*!{suffix}_score.bedgraph.gz \\
            --labels='!{labels_replicates}' \\
            --blacklist='!{blacklist}' |
            gzip -9 > 'insulation_scores_replicates.bedgraph.gz'

        aggregate_insulation_scores.py \\
            *{WT,T1,C1}_*merged!{suffix}_score.bedgraph.gz \\
            --labels='!{labels_conditions}' \\
            --blacklist='!{blacklist}' |
            gzip -9 > 'insulation_scores_conditions.bedgraph.gz'
        '''
}


process generate_tad_report {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    label 'process_short'

    input:
        tuple val(normalization),
              val(resolution),
              path(domains)

        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        generate_tad_report.py \\
            hg38_???_*!{suffix}_domains.bed.gz \\
            --output-prefix=report_replicates \\
            --labels='!{labels_replicates}'

        generate_tad_report.py \\
            *{WT,T1,C1}_*merged!{suffix}_domains.bed.gz \\
            --output-prefix=report_conditions \\
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
        tuple val(normalization),
              val(resolution),
              path(insulation),
              path(blacklist)

        path tads

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        generate_insulation_report.py \\
            insulation_scores_replicates.bedgraph.gz \\
            *WT_*merged!{suffix}_domains.bed.gz \\
            --output-prefix=report_replicates_insulation \\
            --blacklist '!{blacklist}'

        generate_insulation_report.py \\
            insulation_scores_conditions.bedgraph.gz \\
            *WT_*merged!{suffix}_domains.bed.gz \\
            --output-prefix=report_conditions_insulation \\
            --blacklist '!{blacklist}'

        mkdir plots
        mv *.svg *.png plots/
        '''
}

process generate_tad_interaction_scatter {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    input:
        tuple val(normalization),
              val(resolution),
              path(coolers),
              path(tads)

        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        generate_tad_interaction_scatter.py \\
            hg38_???_*!{suffix}.cool \\
            *WT_*merged!{suffix}_domains.bed.gz \\
            --output-prefix=tad_interactions_scatter_replicates \\
            --labels='!{labels_replicates}'

        generate_tad_interaction_scatter.py \\
            *{WT,T1,C1}_merged!{suffix}.cool \\
            *WT_*merged!{suffix}_domains.bed.gz \\
            --output-prefix=tad_interactions_scatter_conditions \\
            --labels='!{labels_conditions}'

        mkdir plots
        mv *.svg *.png plots/
        '''
}

process generate_tad_overlap_report {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "${normalization}/${resolution}/${it}" }

    tag "${normalization}_${resolution}"

    input:
        tuple val(normalization),
              val(resolution),
              path(domains),
              path(blacklist)

        val labels_replicates
        val labels_conditions

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg

    shell:
        suffix="_${normalization}_${resolution}"
        '''
        generate_tad_overlap_report.py \\
            *{WT,T1,C1}_merged!{suffix}_domains.bed.gz \\
            --blacklist '!{blacklist}' \\
            --output-prefix=tad_overlap_conditions \\
            --title='!{normalization}_!{resolution}' \\
            --labels='!{labels_conditions}'

        mkdir plots
        mv *.svg *.png plots/
        '''

}
