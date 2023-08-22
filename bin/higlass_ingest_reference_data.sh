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

sudo -E higlass-manage ingest \
  --filetype chromsizes-tsv \
  --project-name "Chromosomes" \
  --datatype chromsizes \
  --assembly hg38 \
  --name hg38 \
  --uid "$(uuidgen "hg38.filtered.chrom.sizes")" \
  --no-upload \
  --hg-name "$name" \
  input/hg38/hg38.filtered.chrom.sizes
