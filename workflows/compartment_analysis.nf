#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    def mcools = file(params.mcools, checkIfExists: true).sort()
    def resolutions = params.resolutions.sort()
    def dchic_cis_mem = [resolutions, params.dchic_cis_gbs_per_cpu].transpose()

    generate_blacklist(
        file(params.assembly_gaps, checkIfExists: true),
        file(params.cytoband, checkIfExists: true)
    )

    extract_chrom_sizes_from_cooler(
        mcools.first(),
        resolutions.first()
    )

    generate_blacklist.out.bed.set { blacklist }
    extract_chrom_sizes_from_cooler.out.chrom_sizes.set { chrom_sizes }

    filter_fna(
        file(params.ref_genome_fa),
        chrom_sizes
    )

    filter_gtf(
        file(params.ref_gene_gtf),
        chrom_sizes
    )

    filter_fna.out.fna.set { reference_fna }
    filter_gtf.out.gtf.set { refgene_gtf }

    Channel.of(resolutions)
        .flatten()
        .combine(mcools)
        .set { coolers }

    preproc_coolers_for_dchic(
        coolers,
        chrom_sizes,
        blacklist
    )

    preproc_coolers_for_dchic.out.bins
        .groupTuple(sort: true, size: mcools.size())
        .set { bins }

    preproc_coolers_for_dchic.out.biases
        .groupTuple(sort: true, size: mcools.size())
        .set { biases }

    preproc_coolers_for_dchic.out.matrix
        .groupTuple(sort: true, size: mcools.size())
        .set { matrices }

    stage_dchic_inputs(
        biases,
        file(params.dchic_sample_file),
        params.ref_genome_name,
        chrom_sizes,
        reference_fna,
        refgene_gtf
    )

    stage_dchic_inputs.out.tsv
        .join(matrices, failOnDuplicate: true, failOnMismatch: true)
        .join(bins, failOnDuplicate: true, failOnMismatch: true)
        .join(stage_dchic_inputs.out.ref_genome_name, failOnDuplicate: true, failOnMismatch: true)
        .join(stage_dchic_inputs.out.tar, failOnDuplicate: true, failOnMismatch: true)
        .join(Channel.of(dchic_cis_mem).flatten().collate(2), failOnDuplicate: true, failOnMismatch: true)
        .combine([params.dchic_seed])
        .set{ dchic_input_ch }

    run_dchic_cis(dchic_input_ch) |
        run_dchic_select |
        run_dchic_analyze |
        run_dchic_fithic |
        run_dchic_dloop |
        run_dchic_subcomp |
        // run_dchic_enrich |
        run_dchic_viz |
        // 0: resolution
        // 4: input_tar
        // 7: output_prefix
        map { tuple(it[0], it[4], it[7]) } |
        postprocess_dchic_output

    generate_subcompartment_transition_report(
        postprocess_dchic_output.out.subcompartments
    )
    plot_subcompartment_coverage(
        postprocess_dchic_output.out.subcompartments
    )
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

        cat <(zcat '!{assembly_gaps}' | cut -f 2-) \\
            <(zcat '!{cytoband}' | grep 'acen$') |
            grep '^chr[XY0-9]\\+[[:space:]]' |
            cut -f 1-3 |
            sort -k1,1V -k2,2n |
            bedtools merge -i stdin |
            cut -f1-3 > blacklist.bed
        '''
}

process filter_fna {
    label 'process_medium'
    label 'process_short'

    input:
        path fna
        path chrom_sizes

    output:
        path "*fna.gz", emit: fna

    shell:
        outname="${fna.baseName}_filtered.fna.gz"
        '''
        set -o pipefail

        # Filter out chromosomes not listed in the .chrom.sizes file
        zcat '!{fna}' > reference.fa.tmp
        samtools faidx \\
            -r <(cut -f 1 '!{chrom_sizes}') \\
            reference.fa.tmp |
            bgzip -l 9 -@!{task.cpus} > '!{outname}'
        '''
}

process filter_gtf {
    label 'process_very_short'

    input:
        path gtf
        path chrom_sizes

    output:
        path "*gtf.gz", emit: gtf

    shell:
        outname="${gtf.baseName}_filtered.gtf.gz"
        '''
        set -o pipefail

        # Filter out records for chromosomes not listed in the .chrom.sizes file
        awk -F '\\t' 'BEGIN{OFS=FS} {print $1,"0",$2}' '!{chrom_sizes}' |
        bedtools intersect  \\
            -a '!{gtf}'     \\
            -b stdin        \\
            -wa |
            gzip -9 > '!{outname}'
        '''
}

process extract_chrom_sizes_from_cooler {
    label 'process_very_short'

    input:
        path cooler
        val resolution

    output:
        path "*.chrom.sizes", emit: chrom_sizes

    shell:
        outname="${cooler.baseName}_${resolution}.chrom.sizes"
        '''
        set -o pipefail

        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        # Dump .chrom.sizes
        cooler dump -t chroms "$cooler" | grep 'chr[0-9X]\\+' > '!{outname}'
        '''
}

process preproc_coolers_for_dchic {
    label 'process_medium'

    input:
        tuple val(resolution),
              path(cooler)

        path chrom_sizes
        path blacklist

    output:
        tuple val(resolution), path("*_abs.bed.gz"), emit: bins
        tuple val(resolution), path("*.biases.gz"), emit: biases
        tuple val(resolution), path("*.matrix.gz"), emit: matrix

    shell:
        outprefix="${cooler.simpleName}_${resolution}"
        // For some weird reason, adding comments inside the shell block for this process
        // breaks -resume!?!?
        // Comment #1:
        //   Dump the bin table.
        //   4th column contains the absolute bin id (starting from 1)
        //   The 5th column contains 1 if the bin overlaps with one or
        //   more blacklisted regions and 0 otherwise
        // Comment #2:
        //   cooler dump --one-based-ids does not work as intended, so instead we use AWK to offset bin IDs
        '''
        set -o pipefail

        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}'
        else
            cooler='!{cooler}'
        fi

        cooler dump -t bins "$cooler" |
            grep '^chr[0-9X]\\+'      |
            cut -f 1-3                |
            bedtools intersect       \\
                -a stdin             \\
                -b '!{blacklist}'    \\
                -wa -c                |
            awk -F '\\t' 'BEGIN{OFS=FS} {print $1,$2,$3,NR,$4!=0}' |
            pigz -9 -p '!{task.cpus}' > '!{outprefix}_abs.bed.gz'

        dump_bias_vectors_for_fithic.py "$cooler" |
            pigz -9 -p '!{task.cpus}' > '!{outprefix}.biases.gz'

        readarray -t chroms < <(cut -f 1 '!{chrom_sizes}')
        for chrom in "${chroms[@]}"; do
            cooler dump -t pixels "$cooler" -r "$chrom"
        done |
            awk -F '\\t' 'BEGIN{OFS=FS} {print $1+1,$2+1,$3}' |
            pigz -9 -p '!{task.cpus}' > '!{outprefix}.matrix.gz'
        '''
}


process stage_dchic_inputs {
    label 'error_retry'
    label 'process_medium'
    label 'process_short'

    input:
        tuple val(resolution),
              path(biases)

        path sample_file_template

        val ref_genome_name
        path chrom_sizes
        path ref_genome_fa
        path refgene_gtf

    output:
        tuple val(resolution),
              val(ref_genome_name),
              emit: ref_genome_name

        tuple val(resolution),
              path("*.tsv"),
              emit: tsv

        tuple val(resolution),
              path("*.tar.zst"),
              emit: tar

    shell:
        output_prefix="${sample_file_template.simpleName}_${resolution}"
        sample_file="${output_prefix}.tsv"
        ref_folder="${ref_genome_name}_${resolution}_goldenpathData"
        chrom_sizes_=chrom_sizes[0]
        '''
        set -o pipefail

        mkdir biases/ '!{ref_folder}'

        sed 's/{{resolution}}/!{resolution}/g' '!{sample_file_template}' | tee '!{sample_file}'

        for f in *.biases.gz; do
            cp "$f" "biases/$(echo "$f" | sed 's/_!{resolution}\\.biases\\.gz$/.biases.gz/')"
        done

        zcat '!{ref_genome_fa}' > '!{ref_folder}/!{ref_genome_name}.fa'
        cp '!{refgene_gtf}' '!{ref_folder}/!{ref_genome_name}.refGene.gtf.gz'
        cp '!{chrom_sizes_}' '!{ref_folder}/!{ref_genome_name}.chrom.sizes'

        tar -cf - biases/ '!{ref_folder}' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 -o '!{output_prefix}.tar.zst'
        '''
}

process run_dchic_cis {
    label 'error_retry'

    cpus 10
    memory { memory * task.cpus * task.attempt }
    time { 16.hour * task.attempt }


    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              val(ref_genome_name),
              path(input_tar),
              val(memory),
              val(seed)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.cis.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val("${sample_file.simpleName}")

    shell:
        output_prefix="${sample_file.simpleName}"
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS=1

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --cthread '!{task.cpus}' \\
            --seed '!{seed}' \\
            --pcatype cis

        cp -L '!{input_tar}' '!{output_prefix}.cis.tar.zst'
        tar -cf - *_pca/ |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.cis.tar.zst'
        '''
}

process run_dchic_select {
    label 'error_retry'
    label 'process_long'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.select.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --pcatype select

        cp -L '!{input_tar}' '!{output_prefix}.select.tar.zst'
        tar -cf - *_pca/*_pca/*_mat/ |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.select.tar.zst'
        '''
}

process run_dchic_analyze {
    label 'error_retry'
    label 'process_long'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.analyze.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --pcatype analyze

        cp -L '!{input_tar}' '!{output_prefix}.analyze.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.analyze.tar.zst'
        '''
}

process run_dchic_fithic {
    label 'error_retry'
    label 'process_long'
    label 'process_medium'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.fithic.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS=1

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        # If this step fails with an error like:
        #   Error in ids_sample[[i]] : subscript out of bounds
        #   Calls: fithicformat
        #   Execution halted
        # this most likely means the bias vector are not valid (e.g. weights are not centered around 1)
        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --fithicpath="$(which fithic)" \\
            --pythonpath="$(which python3)" \\
            --sthread '!{task.cpus}' \\
            --seed '!{seed}' \\
            --pcatype fithic

        cp -L '!{input_tar}' '!{output_prefix}.fithic.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/fithic_run/' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.fithic.tar.zst'
        '''
}

process run_dchic_dloop {
    label 'error_retry'
    label 'process_long'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.dloop.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --pcatype dloop

        cp -L '!{input_tar}' '!{output_prefix}.dloop.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/fdr_result/' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.dloop.tar.zst'
        '''
}

process run_dchic_subcomp {
    label 'error_retry'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.subcomp.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --pcatype subcomp

        cp -L '!{input_tar}' '!{output_prefix}.subcomp.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/fdr_result/'*.subcompartments.* |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.subcomp.tar.zst'
        '''
}

process run_dchic_enrich {
    label 'error_retry'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.enrich.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --region both \\
            --pcgroup pcQnm \\
            --interaction intra \\
            --compare F \\
            --pcatype enrich

        cp -L '!{input_tar}' '!{output_prefix}.enrich.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/geneEnrichment/' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.enrich.tar.zst'
        '''
}

process run_dchic_viz {
    label 'error_retry'
    label 'process_low'

    input:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path(input_tar),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    output:
        tuple val(resolution),
              path(sample_file),
              path(matrices),
              path(bins),
              path("*.viz.tar.zst"),
              val(ref_genome_name),
              val(seed),
              val(output_prefix)

    shell:
        '''
        set -o pipefail

        mkdir tmp/
        trap 'rm -rf tmp/' EXIT

        export TMPDIR="$PWD/tmp"
        export OPENBLAS_NUM_THREADS='!{task.cpus}'

        zstdcat '!{input_tar}' | tar -xf - --ignore-zeros

        dchicf.r \\
            --file '!{sample_file}' \\
            --dirovwt T \\
            --gfolder *_goldenpathData/ \\
            --genome '!{ref_genome_name}' \\
            --diffdir '!{output_prefix}' \\
            --seed '!{seed}' \\
            --pcatype viz

        cp -L '!{input_tar}' '!{output_prefix}.viz.tar.zst'
        tar -cf - 'DifferentialResult/!{output_prefix}/viz/' |
            zstd -T'!{task.cpus}' --adapt=min=3,max=6 >> '!{output_prefix}.viz.tar.zst'
        '''
}

process postprocess_dchic_output {
    publishDir params.output_dir, mode: 'copy'

    label 'process_high'

    input:
        tuple val(resolution),
              path(input_tar),
              val(output_prefix)

    output:
        tuple val(resolution),
              path("${resolution}/${output_prefix}.dchic.tar.zst"),
              emit: dchic_archive
        tuple val(resolution),
              path("${resolution}/${output_prefix}.pcOri.bedGraph.gz"),
              emit: pc_ori
        tuple val(resolution),
              path("${resolution}/${output_prefix}.Filtered.pcOri.bedGraph.gz"),
              emit: pc_ori_filtered
        tuple val(resolution),
              path("${resolution}/${output_prefix}.pcQnm.bedGraph.gz"),
              emit: pc_qnm
        tuple val(resolution),
              path("${resolution}/${output_prefix}.Filtered.pcQnm.bedGraph.gz"),
              emit: pc_qnm_filtered
        tuple val(resolution),
              path("${resolution}/${output_prefix}.subcompartments.bedGraph.gz"),
              emit: subcompartments
        tuple val(resolution),
              path("${resolution}/${output_prefix}.viz.tar.xz"),
              emit: viz

    shell:
        '''
        set -o pipefail

        mkdir '!{resolution}/' '!{output_prefix}'
        zstdcat '!{input_tar}' |
            tar -C '!{output_prefix}/' -xf - --ignore-zeros \\
                --exclude=biases \\
                --exclude='*_goldenpathData'

        function process_bedgraph {
            bgraph="$1"
            outname="!{output_prefix}$(echo "$bgraph" | sed 's/.*sample_group//')"
            if [[ "$bgraph" == *.subcompartments.* ]]; then
                postprocess_dchic_subcompartments.py "$bgraph" |
                    gzip -9c > "!{resolution}/$outname.gz"
            else
                gzip -9c "$bgraph" > "!{resolution}/$outname.gz"
            fi
        }

        export -f process_bedgraph

        # Compress output files
        srcdir='!{output_prefix}/DifferentialResult/!{output_prefix}/fdr_result'
        printf '%s\\n' "$srcdir/"*sample_group*.bedGraph |
            xargs -P '!{task.cpus}' -I {} bash -c 'process_bedgraph "$@"' _ {}

        tar -C '!{output_prefix}/DifferentialResult' -cf - '!{output_prefix}/viz' |
            xz -T!{task.cpus} -9 --extreme > '!{resolution}/!{output_prefix}.viz.tar.xz'

        tar -cf - '!{output_prefix}' |
            zstd -T!{task.cpus} --adapt=min=13,max=19 -o '!{resolution}/!{output_prefix}.dchic.tar.zst'
        '''
}

process generate_subcompartment_transition_report {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'process_very_short'

    input:
        tuple val(resolution),
              path(bedgraph)

    output:
        val resolution, emit: resolution
        path "${resolution}/plots/*.svg", emit: svg
        path "${resolution}/plots/*.png", emit: png

    shell:
        outprefix="${resolution}/plots/${bedgraph.simpleName}"
        '''
        # Plot alluvials
        generate_compartment_transition_report.py \\
            --plot-type=alluvial \\
            --output-prefix='!{outprefix}_subcompartments' \\
            --base-color="white" \\
            --width=5 \\
            --height=3.5 \\
            --path-to-plotting-script="$(which make_ab_comp_alluvial.r)" \\
            '!{bedgraph}'

        generate_compartment_transition_report.py \\
            --plot-type=alluvial \\
            --aggregate-subcompartments \\
            --output-prefix='!{outprefix}_compartments' \\
            --highlight-color="#b40426ff" \\
            --base-color="#3b4cc0ff" \\
            --path-to-plotting-script="$(which make_ab_comp_alluvial.r)" \\
            '!{bedgraph}'

        # Plot heatmaps
        generate_compartment_transition_report.py \\
            --plot-type=heatmap \\
            --output-prefix='!{outprefix}_subcompartments' \\
            '!{bedgraph}'

        generate_compartment_transition_report.py \\
            --plot-type=heatmap \\
            --aggregate-subcompartments \\
            --output-prefix='!{outprefix}_compartments' \\
            '!{bedgraph}'
        '''
}

process plot_subcompartment_coverage {
    publishDir params.output_dir, mode: 'copy'

    label 'process_low'
    label 'process_very_short'

    input:
        tuple val(resolution),
              path(bedgraph)

    output:
        val resolution, emit: resolution
        path "${resolution}/plots/*.svg", emit: svg
        path "${resolution}/plots/*.png", emit: png
        path "${resolution}/*.tsv", emit: tsv

    shell:
        outprefix="${resolution}/${bedgraph.simpleName}"
        '''
        compare_subcompartment_coverage.py \\
            '!{bedgraph}' \\
            '!{outprefix}_subcompartment_coverage_gw' \\
            --genome-wide

        compare_subcompartment_coverage.py \\
            '!{bedgraph}' \\
            '!{outprefix}_subcompartment_coverage'

        mkdir -p '!{resolution}/plots/'
        mv '!{resolution}/'*.{svg,png} '!{resolution}/plots/'
        '''
}
