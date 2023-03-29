#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    Channel.of("${params.resolutions}".split(','))
           .combine(Channel.fromPath(params.mcools)
                           .filter{ !it.getSimpleName().contains("merged") })
           .set { mcools }

    generate_blacklist(file(params.assembly_gaps),
                       file(params.cytoband))

    preproc_data_for_dchic(mcools,
                           generate_blacklist.out.bed)

    preproc_data_for_dchic.out.resolution
                          .merge(preproc_data_for_dchic.out.matrix)
                          .groupTuple(sort: true)
                          .map { it[1] }
                          .set { matrices }

    preproc_data_for_dchic.out.resolution
                          .merge(preproc_data_for_dchic.out.bins)
                          .groupTuple(sort: true)
                          .map { it[1] }
                          .set { bins }

    preproc_data_for_dchic.out.resolution
                          .merge(preproc_data_for_dchic.out.biases)
                          .groupTuple(sort: true)
                          .map { it[1] }
                          .set { biases }

    preproc_data_for_dchic.out.resolution
                          .merge(preproc_data_for_dchic.out.matrix)
                          .groupTuple(sort: true)
                          .map { it[0] }
                          .set { resolutions }

    run_dchic(matrices,
              bins,
              biases,
              file(params.dchic_sample_file),
              params.ref_genome_name,
              resolutions)

    postprocess_dchic_output(run_dchic.out.resolution,
                             run_dchic.out.pc_ori,
                             run_dchic.out.pc_ori_filtered,
                             run_dchic.out.pc_qnm,
                             run_dchic.out.pc_qnm_filtered,
                             run_dchic.out.subcompartments,
                             run_dchic.out.viz,
                             file(params.dchic_sample_file))

    generate_subcompartment_transition_report(postprocess_dchic_output.out.resolution,
                                              postprocess_dchic_output.out.subcompartments)
    plot_subcompartment_coverage(postprocess_dchic_output.out.resolution,
                                 postprocess_dchic_output.out.subcompartments)
    plot_subcompartment_size_distribution(postprocess_dchic_output.out.resolution,
                                          postprocess_dchic_output.out.subcompartments)
}

process generate_blacklist {
    input:
        path assembly_gaps
        path cytoband

    output:
        path "*.bed", emit: bed

    shell:
        '''
        set -o pipefail

        cat <(zcat '!{assembly_gaps}' | cut -f 2-) \
            <(zcat '!{cytoband}' | grep 'acen$') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 > blacklist.bed
        '''
}

process preproc_data_for_dchic {
    label 'process_low'

    input:
        tuple val(resolution), path(cooler)
        path blacklist

    output:
        path "*.chrom.sizes", emit: chrom_sizes
        path "*.matrix", emit: matrix
        path "*_abs.bed", emit: bins
        path "*.biases.gz", emit: biases
        val resolution, emit: resolution

    shell:
        '''
        set -o pipefail

        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        outprefix="$(basename '!{cooler}' .cool)"
        outprefix="$(basename "$outprefix" .mcool)"

        cooler dump -t chroms "$cooler" | grep 'chr[0-9X]\\+' > "$outprefix.chrom.sizes"

        preprocess.py -input cool                          \
                      -file '!{cooler}'                    \
                      -genomeFile "$outprefix.chrom.sizes" \
                      -res '!{resolution}'                 \
                      -prefix "$outprefix"                 \
                      -removeChr chrY,chrM

        # Dump the bin table.
        # 4th column contains the absolute bin id (starting from 1)
        # The 5th column contains 1 if the bin overlaps with one or more blacklisted regions and 0 otherwise
        cooler dump -t bins "$cooler" |
            grep -v 'chr[0-9X]\\+'    |
            cut -f 1-3                |
            bedtools intersect        \
                -a stdin              \
                -b '!{blacklist}'     \
                -wa -c                |
            awk -F '\\t' 'BEGIN{OFS=FS} {print $1,$2,$3,NR,$4!=0}' |
            tee "${outprefix}_abs.bed" > /dev/null

        cooler dump --na-rep nan -t bins "$cooler" |
            grep -v 'chr[0-9X]\\+'                 |
            awk -F '\\t' 'BEGIN{ OFS=FS } { printf "%s\\t%.0f\\t%s\\n", $1,$2+(!{resolution}/2),$4 }' |
            gzip -9 > "$outprefix.biases.gz"
        '''
}

process run_dchic {
    label 'error_retry'
    label 'process_low'

    memory { task.attempt * 15e9 as Long }  // 15GB

    input:
        path matrices
        path bins
        path biases
        path sample_file

        val ref_genome_name
        val resolution

    output:
        val resolution, emit: resolution
        path "${sample_file.simpleName}_${resolution}.pcOri.bedGraph.gz", emit: pc_ori
        path "${sample_file.simpleName}_${resolution}.Filtered.pcOri.bedGraph.gz", emit: pc_ori_filtered
        path "${sample_file.simpleName}_${resolution}.pcQnm.bedGraph.gz", emit: pc_qnm
        path "${sample_file.simpleName}_${resolution}.Filtered.pcQnm.bedGraph.gz", emit: pc_qnm_filtered
        path "${sample_file.simpleName}_${resolution}.subcompartments.bedGraph.gz", emit: subcompartments
        path "${sample_file.simpleName}_${resolution}.viz.tar.xz", emit: viz

    shell:
        output_prefix="${sample_file.simpleName}_${resolution}"
        '''
        mkdir biases/ tmp/

        TMPDIR="$PWD/tmp"
        export TMPDIR

        sample_file='!{sample_file}'
        sample_file="${sample_file%.tsv}_!{resolution}.tsv"

        sed 's/{{resolution}}/!{resolution}/g' '!{sample_file}' | tee "$sample_file"

        for f in *.biases.gz; do
            ln -sf "$f" "biases/$f"
        done

        dchicf.r --file "$sample_file" --pcatype cis --dirovwt T
        dchicf.r --file "$sample_file" --pcatype select --dirovwt T --genome '!{ref_genome_name}'
        dchicf.r --file "$sample_file" --pcatype analyze --dirovwt T --diffdir '!{output_prefix}'

        # This step fails with an error like: Error in ids_sample[[i]] : subscript out of bounds
        # Skip it for now
        echo \\
        dchicf.r --file "$sample_file" --pcatype fithic --dirovwt T --diffdir '!{output_prefix}' --fithicpath="$(which fithic)" --pythonpath="$(which python3)"

        # Dloop segfaults for resolutions higher than 100kb
        echo \\
        dchicf.r --file "$sample_file" --pcatype dloop --dirovwt T --diffdir '!{output_prefix}'

        dchicf.r --file "$sample_file" --pcatype subcomp --dirovwt T --diffdir '!{output_prefix}'
        dchicf.r --file "$sample_file" --pcatype viz --diffdir '!{output_prefix}' --genome '!{ref_genome_name}'

        # Compress output files
        srcdir='DifferentialResult/!{output_prefix}/fdr_result'
        for f in "$srcdir/"*sample_group*.bedGraph; do
            outname="!{output_prefix}$(echo "$f" | sed 's/.*sample_group//')"
            gzip -9c "$f" > "$outname.gz"
        done

        tar -cf - 'DifferentialResult/!{output_prefix}/viz' |
            xz -T!{task.cpus} -9 > '!{output_prefix}.viz.tar.xz'
        '''
}

process postprocess_dchic_output {
    publishDir params.output_dir, mode: 'copy'

    input:
        val resolution
        path pc_ori
        path pc_ori_filtered
        path pc_qnm
        path pc_qnm_filtered
        path subcompartments
        path viz
        path sample_file

    output:
        val resolution, emit: resolution
        path "${resolution}/${sample_file.simpleName}_${resolution}.pcOri.bedGraph.gz", emit: pc_ori
        path "${resolution}/${sample_file.simpleName}_${resolution}.Filtered.pcOri.bedGraph.gz", emit: pc_ori_filtered
        path "${resolution}/${sample_file.simpleName}_${resolution}.pcQnm.bedGraph.gz", emit: pc_qnm
        path "${resolution}/${sample_file.simpleName}_${resolution}.Filtered.pcQnm.bedGraph.gz", emit: pc_qnm_filtered
        path "${resolution}/${sample_file.simpleName}_${resolution}.subcompartments.bedGraph.gz", emit: subcompartments
        path "${resolution}/${sample_file.simpleName}_${resolution}.viz.tar.xz", emit: viz

    shell:
        '''
        mkdir '!{resolution}/'

        for f in '!{pc_ori}' '!{pc_ori_filtered}' '!{pc_qnm}' '!{pc_qnm_filtered}' '!{viz}'; do
            cp "$f" "!{resolution}/$(basename "$f")"
        done

        outname="!{resolution}/$(basename '!{subcompartments}')"
        postprocess_dchic_subcompartments.py \
            '!{subcompartments}' | gzip -9 > "$outname"
        '''
}

process generate_subcompartment_transition_report {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'process_very_short'

    input:
        val resolution
        path bedgraph

    output:
        val resolution, emit: resolution
        path "${resolution}/plots/*.svg", emit: svg
        path "${resolution}/*.tsv", emit: tsv

    shell:
        outprefix="${resolution}/${bedgraph.simpleName}"
        '''
        generate_compartment_transition_report.py \
            --output-prefix='!{outprefix}_subcompartments' \
            --base-color="white" \
            --width=5 \
            --height=3.5 \
            --path-to-plotting-script="$(which make_ab_comp_alluvial.r)" \
            '!{bedgraph}'

        generate_compartment_transition_report.py \
            --aggregate-subcompartments \
            --output-prefix='!{outprefix}_compartments' \
            --highlight-color="#b40426ff" \
            --base-color="#3b4cc0ff" \
            --path-to-plotting-script="$(which make_ab_comp_alluvial.r)" \
            '!{bedgraph}'

        mkdir -p '!{resolution}/plots/'
        mv '!{resolution}/'*.svg '!{resolution}/plots/'
        '''
}

process plot_subcompartment_coverage {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'process_very_short'

    input:
        val resolution
        path bedgraph

    output:
        val resolution, emit: resolution
        path "${resolution}/plots/*.svg", emit: svg
        path "${resolution}/plots/*.png", emit: png
        path "${resolution}/*.tsv", emit: tsv

    shell:
        outprefix="${resolution}/${bedgraph.simpleName}"
        '''
        compare_subcompartment_coverage.py \
            '!{bedgraph}' \
            '!{outprefix}_subcompartment_coverage_gw' \
            --genome-wide

        compare_subcompartment_coverage.py \
            '!{bedgraph}' \
            '!{outprefix}_subcompartment_coverage'

        mkdir -p '!{resolution}/plots/'
        mv '!{resolution}/'*.{svg,png} '!{resolution}/plots/'
        '''
}

process plot_subcompartment_size_distribution {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'process_very_short'

    input:
        val resolution
        path bedgraph

    output:
        val resolution, emit: resolution
        path "${resolution}/plots/*.svg", emit: svg
        path "${resolution}/plots/*.png", emit: png

    shell:
        outprefix="${resolution}/${bedgraph.simpleName}_size_distribution"
        '''
        compare_subcompartment_size_distribution.py \
            '!{bedgraph}' \
            '!{outprefix}'

        mkdir -p '!{resolution}/plots/'
        mv '!{resolution}/'*.{svg,png} '!{resolution}/plots/'
        '''
}
