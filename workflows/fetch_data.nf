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
    checksums = file(params.checksums)

    // Read download_list file line by line and store the file content in the urls list
    Channel.fromPath(dls)
            .splitText()
            .map {
                 // Lines are formatted like "url\tdest_name\tdest_dir"
                 toks = it.trim().split('\t')
                 assert toks.size() == 3

                 url = file(toks[0])

                 base_name = normalize_path_name(url.getName().toString())
                 dest_name = toks[1]
                 dest_dir = toks[2]

                 tuple(url, base_name, dest_dir, dest_name)
                 }
            .set { urls }

    validate_files(urls, checksums)

    rename_and_compress_files(validate_files.out.file, validate_files.out.metadata)
}

process validate_files {
    label 'process_short'

    input:
        tuple path(url), val(base_name), val(dest_dir), val(dest_name)
        path checksums

    output:
        tuple val(base_name), val(dest_dir), val(dest_name), emit: metadata
        path "*.ok", emit: file

    shell:
        src="${url.fileName}"
        dest=normalize_path_name(src)
        '''
        sha256sum -c '!{checksums}' --ignore-missing
        mv '!{src}' '!{dest}.ok'
        '''
}

process rename_and_compress_files {
    publishDir "${params.data_dir}", mode: 'copy',
                                       saveAs: { "${dest_dir}/${dest_name}" }
    label 'process_short'

    input:
        path src
        tuple val(base_name), val(dest_dir), val(dest_name)

    output:
        path "${dest_dir}/${dest_name}.new", emit: file

    shell:
        '''
        if='!{src}'
        of='!{dest_dir}/!{dest_name}.new'

        mkdir -p '!{dest_dir}'

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