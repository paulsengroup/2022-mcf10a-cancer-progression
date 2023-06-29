#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


def collect_files(prefix, sample_id, suffix, type = "file") {
    def files = file("${prefix}/${sample_id}*${suffix}",
                     type: type,
                     checkIfExists: true)

    // Example: hg38_001_MCF10A_WT_REP1 -> hg38_MCF10A_WT
    def condition_id = sample_id.replaceAll(/(.*)_\d{3}_(.*)_REP\d/, '$1_$2')

    return tuple(sample_id, condition_id, files)
}


workflow {
    Channel.fromPath(params.nfcore_samplesheet)
        .splitCsv(sep: ",", header: true)
        .map { it.sample }
        .unique()
        .set { sample_ids }

    input_dir = file(params.input_dir, type: 'dir', checkIfExists: true)

    sample_ids
        .map { collect_files("${input_dir}/pairs/", it, '.allValidPairs.xz') }
        .set { input_valid_pairs }

    Channel.fromPath("${input_dir}/stats/*.tar.xz", checkIfExists: true, type: 'file')
        .set { stat_tars }

    cooler_cload(input_valid_pairs,
                 file(params.chrom_sizes),
                 params.assembly_name,
                 params.resolutions.first())

    cooler_cload.out.cool.set { coolers_by_sample }

    // Group coolers by condition
    coolers_by_sample
        .groupTuple(by: 1, size: 2)
        // [condition_id, [sample_ids], [coolers]]
        .set { coolers_by_condition }

    cooler_merge(coolers_by_condition)

    Channel.empty()
        .mix(coolers_by_sample.map { tuple(it[0], it[2]) },      // [sample_id, cooler]
             cooler_merge.out.cool.map { tuple(it[0], it[2]) })  // [condition_id, cooler]
        .set { coolers }

    cooler_to_hic(coolers,
                  params.resolutions.join(" "))

    cooler_zoomify(coolers)
    cooler_balance(cooler_zoomify.out.mcool,
                   file(params.blacklist, checkIfExists: true))

    // NOTE label and sample mappings is guaranteed to be correct
    // regardless of the order in which items are emitted by compress_stats()
    // because plot_stats is globbing file names, thus sorting files by their
    // sample number
    plot_stats(stat_tars.collect(),
               params.plot_pretty_labels)
}


process cooler_cload {
    tag "$sample"

    input:
        tuple val(sample),
              val(condition),
              path(pairs)

        val chrom_sizes
        val assembly
        val resolution

    output:
        tuple val(sample),
              val(condition),
              path("*.cool"),
        emit: cool

    shell:
        outname = "${pairs.simpleName}.cool"
        '''
        chrom_sizes="$(mktemp chrom.sizes.XXXXXX)"

        grep '^chr[XY0-9]\\+[[:space:]]' '!{chrom_sizes}' > "$chrom_sizes"

        xz -dcf '!{pairs}' |
        cooler cload pairs \\
            "$chrom_sizes:!{resolution}" \\
            -                            \\
            '!{outname}'                 \\
            --chrom1 2 --pos1 3          \\
            --chrom2 5 --pos2 6          \\
            --assembly='!{assembly}'
        '''
}

process cooler_merge {
    tag "$condition"

    input:
        tuple val(samples),
              val(condition),
              path(coolers)

    output:
        tuple val("${condition}_merged"),
              val(samples),
              path("*.cool"),
        emit: cool

    shell:
        outname = "${condition}_merged.cool"
        '''
        cooler merge '!{outname}' *.cool
        '''
}

process cooler_to_hic {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "hic/${label}.hic" }

    label 'process_medium'
    tag "$label"

    input:
        tuple val(label),
              path(cool)

        val resolutions

    output:
        tuple val(label),
              path("*.hic"),
        emit: hic

    shell:
        memory_gb=task.memory.toGiga()
        out="${cool.baseName}.hic"
        '''
        trap 'rm -rf tmp/' EXIT

        hic_tools_jar=(/usr/local/share/java/hic_tools/hic_tools*.jar)
        juicer_tools_jar=(/usr/local/share/java/juicer_tools/juicer_tools*.jar)

        cool2hic-ng \\
            --hic-tools-jar "${hic_tools_jar[*]}"       \\
            --juicer-tools-jar "${juicer_tools_jar[*]}" \\
            --tmpdir=tmp/                               \\
            '!{cool}'                                   \\
            '!{out}'                                    \\
            --nproc '!{task.cpus}'                      \\
            -Xmx '!{memory_gb}g'                        \\
            --resolutions !{resolutions}
        '''
}

process cooler_zoomify {
    label 'process_medium'
    label 'process_long'
    tag "$label"

    input:
        tuple val(label),
              path(cool)

    output:
        tuple val(label),
              path("*.mcool"),
        emit: mcool

    shell:
        '''
        cooler zoomify -p !{task.cpus} \\
                       -r N            \\
                       '!{cool}'
        '''
}

process cooler_balance {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { "mcools/${label}.mcool" }

    label 'process_medium'
    label 'process_long'
    tag "$label"

    input:
        tuple val(label),
              path(mcool)
        path blacklist

    output:
        tuple val(label),
              path("*.mcool.new"),
        emit: mcool

    shell:
        memory_gb=task.memory.toGiga()
        '''
        hic_tools_jar=(/usr/local/share/java/hic_tools/hic_tools*.jar)
        juicer_tools_jar=(/usr/local/share/java/juicer_tools/juicer_tools*.jar)

        cp '!{mcool}' '!{mcool}.new'
        cooler_balance.py \\
            --hic-tools-jar "${hic_tools_jar[*]}"       \\
            --juicer-tools-jar "${juicer_tools_jar[*]}" \\
            --nproc !{task.cpus}                        \\
            --blacklist '!{blacklist}'                  \\
            -Xmx !{memory_gb}g                          \\
            '!{mcool}'                                  \\
            '!{mcool}.new'
        '''
}

process plot_stats {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path archives
        val pretty_labels

    output:
        path "plots/*.png", emit: png
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        '''
        generate_hic_mapping_report.py \\
            *.tar.xz \\
            hic_mapping_report \\
            --sample-labels '!{pretty_labels}'

        mkdir plots
        mv *.png *.svg plots/
        '''
}
