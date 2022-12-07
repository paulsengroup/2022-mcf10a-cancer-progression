#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    call_compartments(Channel.fromPath(params.mcools),
                      file(params.ref_genome),
                      params.resolution)

    preproc_data_for_dchic(Channel.fromPath(params.mcools)
                                  .filter{ !it.getSimpleName().contains("merged") },
                           params.resolution)

    run_dchic(preproc_data_for_dchic.out.matrix.collect(),
              preproc_data_for_dchic.out.bins.collect(),
              preproc_data_for_dchic.out.biases.collect(),
              Channel.fromPath(params.dchic_sample_files),
              params.ref_genome_name,
              params.resolution)

    plot_subcompartment_transitions(run_dchic.out.subcompartments)
    plot_subcompartment_coverage(run_dchic.out.subcompartments)
    plot_subcompartment_size_distribution(run_dchic.out.subcompartments)
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
    publishDir "${params.output_dir}/diff_compartments", mode: 'copy'

    label 'process_low'

    input:
        path matrices
        path bins
        path biases
        path sample_file

        val ref_genome_name
        val resolution

    output:
        path "*.pcOri.bedGraph.gz", emit: pc_ori
        path "*.Filtered.pcOri.bedGraph.gz", emit: pc_ori_filtered
        path "*.pcQnm.bedGraph.gz", emit: pc_qnm
        path "*.Filtered.pcQnm.bedGraph.gz", emit: pc_qnm_filtered
        path "*.subcompartments.bedGraph.gz", emit: subcompartments
        path "*.viz.tar.xz", emit: viz

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

process plot_subcompartment_transitions {
    publishDir "${params.output_dir}/diff_compartments/", mode: 'copy'

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
        '!{params.script_dir}/plot_compartment_transitions.py' \
            --output-prefix='!{outprefix}_subcompartments' \
            --path-to-plotting-script='!{params.script_dir}/make_ab_comp_alluvial.r' \
            '!{bedgraph}'

        '!{params.script_dir}/plot_compartment_transitions.py' \
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
    publishDir "${params.output_dir}/diff_compartments/", mode: 'copy'

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
    publishDir "${params.output_dir}/diff_compartments/", mode: 'copy'

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
