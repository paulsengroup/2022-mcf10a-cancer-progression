#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2


workflow {
    // Read download_list file line by line and store the file content in downloads
    Channel.fromPath(file(params.download_list))
           .splitCsv(sep: "\t", header: true)
           .map { row -> out_name = file(row.dest).getName()
                         // Keep the path relative
                         out_dir = row.dest.replaceFirst(/${out_name}$/, "")
                         tuple(row.url,
                               row.sha256,
                               out_dir,
                               out_name,
                               row.compress)
                }
           .set { downloads }

    process_files(downloads)

    // Collect checksums and sort them by file path
    process_files.out.checksum
                 .collectFile(name: "checksums.sha256",
                              storeDir: params.data_dir,
                              newLine: true,
                              sort: { it.split("  ")[1] })

}

process process_files {
    publishDir "${params.data_dir}", mode: 'copy',
                                     saveAs: { "${out_dir}/${out_name}" }
    label 'process_short'

    input:
        tuple val(url), val(checksum), val(out_dir), val(out_name), val(compress)

    output:
        path "*.tmp", emit: file
        stdout emit: checksum

    shell:
        dest="${out_name}.tmp"
        final_dest="${out_dir}/${out_name}".replaceAll(/\/\//, "/")  // relpace //
        '''
        set -o pipefail

        hash='!{checksum}'
        tmp_file="$(mktemp -t '!{dest}.XXXXXXXXXX')"
        printf '%s  %s' "$hash" "$tmp_file" > checksum.sha256

        trap 'rm -f "$tmp_file"' EXIT

        # I am using curl instead of letting Nextflow handle remote files
        # because the latter sometime truncates files
        curl -L '!{url}' -o "$tmp_file"
        if ! sha256sum --quiet -c checksum.sha256 1>&2 ; then
            2>& echo "Checksum failed for file \\"!{dest}\\" (\\"!{url}\\")"
            exit 1
        fi

        # Rename and compress files when appropriate
        if [[ '!{compress}' == TRUE ]]; then
            pigz -p !{task.cpus} -9c "$tmp_file" > '!{dest}'
            hash="$(sha256sum '!{dest}' | cut -d' ' -f 1)"
        else
            mv "$tmp_file" '!{dest}'
        fi

        printf '%s  %s' "$hash" '!{final_dest}'
        '''
}
