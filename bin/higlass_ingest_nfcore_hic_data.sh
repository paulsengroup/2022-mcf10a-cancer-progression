#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"

data_dir="$git_root/data/"
higlass_dir="$data_dir/higlass/"

# rm -rf "$higlass_dir"

"$(dirname "$0")/higlass_start_instance.sh"

function uuidgen {
    "$(dirname "$0")/higlass_uuid_generator.py" "$1"
}

export uuidgen

for mcool in "$data_dir/output/nfcore_hic/mcools/GRCh38_0"??_*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "hic_by_sample" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$(basename "$mcool")")" \
        --hg-name "$name" \
        --no-upload \
        --assembly hg38 \
        "${mcool#"$data_dir/"}"
done

for mcool in "$data_dir/output/nfcore_hic/mcools/"*_merged.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "hic_by_condition" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$(basename "$mcool")")" \
        --hg-name "$name" \
        --no-upload \
        --assembly hg38 \
        "${mcool#"$data_dir/"}"
done
