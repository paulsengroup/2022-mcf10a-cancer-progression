#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
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

for beddb in "$data_dir/output/tad_analysis/clodius/GRCh38_0"??_*.beddb; do
    beddb_name="$(basename "$beddb" _domains.bed.beddb)"
    sudo -E higlass-manage ingest \
        --project-name "tad_by_sample" \
        --name "$beddb_name" \
        --uid "$(uuidgen "$(basename "$beddb")")" \
        --hg-name "$name" \
        --no-upload \
        --assembly hg38 \
        "${beddb#"$data_dir/"}"
done

for beddb in "$data_dir/output/tad_analysis/clodius/GRCh38_"*merged*.beddb; do
    beddb_name="$(basename "$beddb" _domains.bed.beddb)"
    sudo -E higlass-manage ingest \
        --project-name "tad_by_condition" \
        --name "$beddb_name" \
        --uid "$(uuidgen "$(basename "$beddb")")" \
        --hg-name "$name" \
        --no-upload \
        --assembly hg38 \
        "${beddb#"$data_dir/"}"
done
