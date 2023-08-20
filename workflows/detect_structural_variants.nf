#!/usr/bin/env nextflow
// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

def collect_files(prefix, sample_id, suffix, type = "file") {
    def files = file("${prefix}/${sample_id}*${suffix}",
                     type: type,
                     checkIfExists: true)

    // Example: hg38_001_MCF10A_WT_REP1 -> hg38_MCF10A_WT
    def condition_id = sample_id.replaceAll(/(.*)_\d{3}_(.*)_REP\d/, '$1_$2')

    return tuple(condition_id, files)
}


workflow {
    Channel.fromPath(params.mcools, checkIfExists: true)
        .set { mcools }

    Channel.fromPath(params.nfcore_hic_samplesheet)
        .splitCsv(sep: ",", header: true)
        .map { it.sample }
        .unique()
        .set { sample_ids }

    aln_dir = file(params.alignment_dir, type: 'dir', checkIfExists: true)

    sample_ids
        .map { collect_files(aln_dir, it, '.bwt2pairs.cram') }
        .set { alignments }

    extract_chrom_sizes_from_cooler(
        mcools.first(),
        params.hictrans_resolution
    )

    extract_chrom_sizes_from_cooler.out.chrom_sizes
        .set { chrom_sizes }

    generate_chromosome_pairs(
        chrom_sizes
    )

    generate_chromosome_pairs.out.pairs
        .splitCsv(sep: '\t', header: true)
        .set { chrom_pairs }

    preproc_coolers_for_hictrans(
        mcools.combine([params.hictrans_resolution])
            .combine(chrom_pairs)
            .map { tuple(it[0], it[1], it[2].chrom1, it[2].chrom2) }
    )

    run_hictrans(
        preproc_coolers_for_hictrans.out.data,
        chrom_sizes
    )

    digest_genome_for_hint(
        file(params.reference_genome, checkIfExists: true),
        params.ref_genome_name,
        params.restriction_enzymes,
        params.restriction_enzymes_alias
    )

    run_hint_cnv(
        mcools,
        file(params.hint_refzip, checkIfExists: true),
        digest_genome_for_hint.out.txt,
        params.hint_resolution,
        params.ref_genome_name
    )

    run_hint_tl(
        mcools,
        file(params.hint_refzip, checkIfExists: true),
        file(params.hint_backdirzip, checkIfExists: true),
        digest_genome_for_hint.out.txt,
        params.ref_genome_name
    )

    filter_mappings(
        alignments,
        params.hic_breakfinder_quality_score
    )

    merge_bams(
        filter_mappings.out.bam.groupTuple()
    )

    run_hic_breakfinder(
        merge_bams.out.bam,
        file(params.hic_breakfinder_expected_intra),
        file(params.hic_breakfinder_expected_inter)
    )
}


process generate_chromosome_pairs {
    label 'process_very_short'

    input:
        path chrom_sizes

    output:
        path "*.tsv", emit: pairs

    shell:
        '''
        #!/usr/bin/env python3

        import pandas as pd

        df = pd.read_table("!{chrom_sizes}", names=["chrom", "length"])

        chroms = df["chrom"].tolist()

        with open("chromosome_pairs.tsv", "w") as f:
            print("chrom1\\tchrom2", file=f)
            for i, chrom1 in enumerate(chroms):
                for chrom2 in chroms[i:]:
                    print(f"{chrom1}\\t{chrom2}", file=f)
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

process preproc_coolers_for_hictrans {
    label 'process_medium'
    tag "${cooler.simpleName}_${chrom1}_${chrom2}"

    input:
        tuple path(cooler),
              val(resolution),
              val(chrom1),
              val(chrom2)

    output:
        tuple val(chrom1),
              val(chrom2),
              path("*_abs.bed.gz"),
              path("*.matrix.gz"), emit: data

    shell:
        outprefix="${cooler.simpleName}"
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
            awk -F '\\t' 'BEGIN{OFS=FS} {print $1,$2,$3,NR}' |
            pigz -9 -p '!{task.cpus}' > '!{outprefix}_abs.bed.gz'

        cooler dump -t pixels "$cooler" -r '!{chrom1}' -r2 '!{chrom2}' |
            awk -F '\\t' 'BEGIN{OFS=FS} {print $1+1,$2+1,$3}' |
            pigz -9 -p '!{task.cpus}' > '!{outprefix}.matrix.gz'
        '''
}


process run_hictrans {
    publishDir "${params.output_dir}/hictrans", mode: 'copy'
    tag "${matrix.simpleName}_${chrom1}_${chrom2}"

    label 'process_long'

    input:
        tuple val(chrom1),
              val(chrom2),
              path(bins),
              path(matrix)

        path chrom_sizes

    output:
        path "*.tar.gz", emit: tar

    shell:
        prefix="${matrix.simpleName}"
        outdir="${prefix}_${chrom1}_${chrom2}"
        '''
        hictrans \\
            --mat '!{matrix}' \\
            --bed '!{bins}' \\
            --chrA '!{chrom1}' \\
            --chrB '!{chrom2}' \\
            --resolutions 2,3,4,5,6,8,10 \\
            --covq 0.1 \\
            --chromsize '!{chrom_sizes}' \\
            --prefix '!{prefix}'

        mkdir -p '!{outdir}'

        if compgen -G "./*/*/MultiResolution_supported_Translocations/*.txt" > /dev/null; then
            cp ./*/*/MultiResolution_supported_Translocations/*.txt '!{outdir}/'
        fi

        if compgen -G "./*/*/Translocations/Details/*.txt" > /dev/null; then
            cp ./*/*/Translocations/Details/*.txt '!{outdir}/'
        fi

        tar -czf '!{outdir}.tar.gz' '!{outdir}'
        '''
}

process digest_genome_for_hint {

    input:
        path fna
        val assembly
        val enzymes
        val enzyme_alias

    output:
        path "*.txt", emit: txt

    shell:
        outputname="${assembly}_${enzyme_alias}_enzymeSites.txt"
        '''
        compute_restriction_sites_for_hint.py \\
            '!{fna}' \\
            !{enzymes} > '!{outputname}'
        '''
}

process run_hint_cnv {
    publishDir "${params.output_dir}/hint_cnv", mode: 'copy'

    tag "${cooler.simpleName}"

    input:
        path cooler
        path refzip
        path restriction_sites
        val resolution
        val ref_genome_name

    output:
        path "*.tar.gz", emit: tar

    shell:
        outprefix="${cooler.simpleName}"
        '''
        if [ '!{resolution}' -ne 0 ]; then
            cooler='!{cooler}::/resolutions/!{resolution}000'
        else
            cooler='!{cooler}'
        fi

        unzip '!{refzip}'
        cp '!{restriction_sites}' hg38/

        hint cnv \\
            --matrixfile "$cooler" \\
            --format cooler \\
            --resolution '!{resolution}' \\
            --refdir '!{ref_genome_name}' \\
            --genome '!{ref_genome_name}' \\
            --name '!{outprefix}' \\
            --bicseq /opt/BICseq2-seg \\
            --enzyme 'arima'

        mv HiNTcnv_OUTPUT '!{outprefix}'
        tar -czf '!{outprefix}.tar.gz' '!{outprefix}'
        '''
}


process run_hint_tl {
    publishDir "${params.output_dir}/hint_tl", mode: 'copy'

    tag "${mcool.simpleName}"

    label 'process_high'

    input:
        path mcool
        path refzip
        path backdirzip
        path restriction_sites
        val ref_genome_name

    output:
        path "*.tar.gz", emit: tar

    shell:
        outprefix="${mcool.simpleName}"
        '''
        unzip '!{refzip}'
        unzip '!{backdirzip}'
        cp '!{restriction_sites}' hg38/

        # Setting --name breaks things
        hint tl \\
            --matrixfile '!{mcool}::/resolutions/1000000,!{mcool}::/resolutions/100000' \\
            --format cooler \\
            --refdir '!{ref_genome_name}' \\
            --backdir '!{ref_genome_name}' \\
            --genome '!{ref_genome_name}' \\
            --cutoff 0.05 \\
            --enzyme 'arima' \\
            --ppath "$(which pairix)" \\
            --threads '!{task.cpus}'

        mv HiNTtransl_OUTPUT '!{outprefix}'
        tar -czf '!{outprefix}.tar.gz' '!{outprefix}'
        '''
}

process filter_mappings {
    label 'process_medium'

    tag "${aln.simpleName}"

    input:
        tuple val(condition),
              path(aln)

        val score

    output:
        tuple val(condition),
              path("*.bam"),
        emit: bam

    shell:
        outname="${aln.baseName}.filtered.bam"
        '''
        mkdir -p tmp
        samtools view -bSq '!{score}' -@'!{task.cpus}' '!{aln}' |
           samtools sort -T tmp/ -O bam -@'!{task.cpus}' -o '!{outname}'
        '''
}

process merge_bams {
    label 'process_medium'

    tag "${condition}"

    input:
        tuple val(condition),
              path(bams)

    output:
        tuple val(condition),
              path("*.bam"),
        emit: bam

    shell:
        outname="${condition}.filtered.bam"
        '''
        samtools merge  -@'!{task.cpus}' -O bam -o '!{outname}' !{bams}
        '''
}


process run_hic_breakfinder {
    publishDir "${params.output_dir}/hic_breakfinder", mode: 'copy'

    tag "${condition}"

    input:
        tuple val(condition),
              path(bam)

        path expected_intra
        path expected_inter

    output:
        path "*breaks*.txt", emit: txt

    shell:
        outprefix=condition
        '''
        hic_breakfinder \\
            --bam-file '!{bam}' \\
            --exp-file-inter '!{expected_inter}' \\
            --exp-file-intra '!{expected_intra}' \\
            --min-1kb \\
            --name '!{outprefix}'
        '''
}
