#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    preproc_data_for_dchic(Channel.fromPath(params.mcools)
                                  .filter{ !it.getSimpleName().contains("merged") },
                           params.resolution)

    run_dchic(preproc_data_for_dchic.out.matrix.collect(),
              preproc_data_for_dchic.out.bins.collect(),
              preproc_data_for_dchic.out.biases.collect(),
              file(params.dchic_sample_file),
              params.ref_genome_name,
              params.resolution)

    postprocess_dchic_output(run_dchic.out.pc_ori,
                             run_dchic.out.pc_ori_filtered,
                             run_dchic.out.pc_qnm,
                             run_dchic.out.pc_qnm_filtered,
                             run_dchic.out.subcompartments,
                             run_dchic.out.viz,
                             file(params.dchic_sample_file))

    generate_subcompartment_transition_report(postprocess_dchic_output.out.subcompartments)
    plot_subcompartment_coverage(postprocess_dchic_output.out.subcompartments)
    plot_subcompartment_size_distribution(postprocess_dchic_output.out.subcompartments)
}

process preproc_data_for_dchic {
    label 'process_low'

    input:
        path cooler
        val resolution

    output:
        path "*.chrom.sizes", emit: chrom_sizes
        path "*.matrix", emit: matrix
        path "*_abs.bed", emit: bins
        path "*.biases.gz", emit: biases

    shell:
        '''
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

        cooler dump --na-rep nan -t bins "$cooler" |
            grep -v 'chr[0-9X]\\+'                 |
            awk -F '\\t' 'BEGIN{ OFS=FS } { printf "%s\\t%.0f\\t%s\\n", $1,$2+(!{resolution}/2),$4 }' |
            gzip -9 > "$outprefix.biases.gz"
        '''
}

process run_dchic {
    label 'process_low'

    input:
        path matrices
        path bins
        path biases
        path sample_file

        val ref_genome_name
        val resolution

    output:
        path "${sample_file.simpleName}.pcOri.bedGraph.gz", emit: pc_ori
        path "${sample_file.simpleName}.Filtered.pcOri.bedGraph.gz", emit: pc_ori_filtered
        path "${sample_file.simpleName}.pcQnm.bedGraph.gz", emit: pc_qnm
        path "${sample_file.simpleName}.Filtered.pcQnm.bedGraph.gz", emit: pc_qnm_filtered
        path "${sample_file.simpleName}.subcompartments.bedGraph.gz", emit: subcompartments
        path "${sample_file.simpleName}.viz.tar.xz", emit: viz

    shell:
        output_prefix="${sample_file.simpleName}"
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
    publishDir params.output_dir, mode: 'copy',
                                  saveAs: { file(it).getName() }

    input:
        path pc_ori
        path pc_ori_filtered
        path pc_qnm
        path pc_qnm_filtered
        path subcompartments
        path viz
        path sample_file

    output:
        path "out/${sample_file.simpleName}.pcOri.bedGraph.gz", emit: pc_ori
        path "out/${sample_file.simpleName}.Filtered.pcOri.bedGraph.gz", emit: pc_ori_filtered
        path "out/${sample_file.simpleName}.pcQnm.bedGraph.gz", emit: pc_qnm
        path "out/${sample_file.simpleName}.Filtered.pcQnm.bedGraph.gz", emit: pc_qnm_filtered
        path "out/${sample_file.simpleName}.subcompartments.bedGraph.gz", emit: subcompartments
        path "out/${sample_file.simpleName}.viz.tar.xz", emit: viz

    shell:
        '''
        mkdir out/

        for f in '!{pc_ori}' '!{pc_ori_filtered}' '!{pc_qnm}' '!{pc_qnm_filtered}' '!{viz}'; do
            cp "$f" "out/$(basename "$f")"
        done

        outname="out/$(basename '!{subcompartments}')"
        '!{params.script_dir}/postprocess_dchic_subcompartments.py' \
            '!{subcompartments}' | gzip -9 > "$outname"
        '''
}

process generate_subcompartment_transition_report {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'very_short'

    input:
        path bedgraph

    output:
        path "plots/*.svg", emit: svg
        path "*.tsv", emit: tsv

    shell:
        outprefix="${bedgraph.simpleName}"
        '''
        '!{params.script_dir}/generate_compartment_transition_report.py' \
            --output-prefix='!{outprefix}_subcompartments' \
            --base-color="white" \
            --width=5 \
            --height=3.5 \
            --path-to-plotting-script='!{params.script_dir}/make_ab_comp_alluvial.r' \
            '!{bedgraph}'

        '!{params.script_dir}/generate_compartment_transition_report.py' \
            --aggregate-subcompartments \
            --output-prefix='!{outprefix}_compartments' \
            --highlight-color="#b40426ff" \
            --base-color="#3b4cc0ff" \
            --path-to-plotting-script='!{params.script_dir}/make_ab_comp_alluvial.r' \
            '!{bedgraph}'

        mkdir plots/
        mv *.svg plots/
        '''
}

process plot_subcompartment_coverage {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'very_short'

    input:
        path bedgraph

    output:
        path "plots/*.svg", emit: svg
        path "plots/*.png", emit: png
        path "*.tsv", emit: tsv

    shell:
        outprefix="${bedgraph.simpleName}"
        '''
        '!{params.script_dir}/compare_subcompartment_coverage.py' \
            '!{bedgraph}' \
            '!{outprefix}_subcompartment_coverage_gw' \
            --genome-wide

        '!{params.script_dir}/compare_subcompartment_coverage.py' \
            '!{bedgraph}' \
            '!{outprefix}_subcompartment_coverage'

        mkdir plots/
        mv *.svg  *.png plots/
        '''
}

process plot_subcompartment_size_distribution {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'very_short'

    input:
        path bedgraph

    output:
        path "plots/*.svg", emit: svg
        path "plots/*.png", emit: png

    shell:
        outprefix="${bedgraph.simpleName}_size_distribution"
        '''
        '!{params.script_dir}/compare_subcompartment_size_distribution.py' \
            '!{bedgraph}' \
            '!{outprefix}'

        mkdir plots/
        mv *.svg  *.png plots/
        '''
}
