#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"

data_dir="$git_root/data/"
higlass_dir="$data_dir/higlass/"

mkdir -p "$higlass_dir"

sudo -E higlass-manage start \
    -d "$higlass_dir" \
    -m "$data_dir" \
    -n "$name" \
    --use-redis \
    --no-public-data \
    -p 9000

sudo chown -R "$USER:$USER" "$higlass_dir"
find "$higlass_dir" -type d -exec chmod 755 {} +
find "$higlass_dir" -type f ! -executable -exec chmod 644 {} +
