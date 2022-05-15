#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

// For some reason importing this function from utils.nfm causes an error like:
//   Missing process or function with name 'normalize_path_name' -- Did you mean 'normalize_path_name' instead?

// Remove problematic characters from file names.
def normalize_path_name(old_path) {
   old_path.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

workflow {
    dls = file(params.download_list)
    chksums = file(params.checksums)

    // Read download_list file line by line and store the file content in the urls list
    Channel.fromPath(dls)
            .splitText()
            .map {
                 // Lines are formatted like "url\tshort_name"
                 toks = it.trim().split('\t')
                 file(toks[0])
                 }
            .set { files }

    validate_files(files, chksums)

    def fn_mappings = [:]
    file(dls)
        .eachLine { line ->
            // Lines are formatted like "url\tshort_name"
            toks = line.trim().split('\t')
            old_name = normalize_path_name(file(toks[0]).getName().toString())
            new_name = file(toks[1]).getName()
            // println "${old_name} -> ${new_name}"
            fn_mappings[old_name] = new_name
    }

    rename_and_compress_files(validate_files.out.files, fn_mappings)
}

process validate_files {
    publishDir "${params.download_dir}", mode: 'copy'

    label 'process_short'

    input:
        path download
        path checksums

    output:
        path "*.ok", emit: files

    shell:
        input_name = "${download.fileName}"
        out_name = normalize_path_name(input_name)
        '''
        sha256sum -c '!{checksums}' --ignore-missing
        mv '!{input_name}' '!{out_name}.ok'
        '''
}

process rename_and_compress_files {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { fname ->
                                                 old_name = file(fname).getBaseName() // Trim .new
                                                 assert name_mappings.get(old_name) != null
                                                 name_mappings[old_name]
                                               }
    label 'process_short'

    input:
        path old_name
        val name_mappings

    output:
        path "*.new"
        val 1, emit: flag

    shell:
        '''
        if="!{old_name}"
        of="${if%.ok}.new"

        # Rename files and compress them
        if [[ $if == *.gz.ok    ||
              $if == *.hic.ok   ||
              $if == *.mcool.ok ||
              $if == *.bigWig.ok ]]; then
            cp "$if" "$of"
        else
            pigz -9c "$if" > "$of"
        fi
        '''
}