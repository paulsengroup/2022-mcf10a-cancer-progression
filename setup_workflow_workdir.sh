#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u
# set -x


if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 repository-root working-dir"
  2>&1 echo "Example: $0 . .nextflow-fetch_data-wd/"
  exit 1
fi


src_dir="$1"
dest_dir="$2"

if [ ! -d "$src_dir" ]; then
  2>&1 echo "Repository root does not exist: '$src_dir'"
  exit 1
fi


mkdir -p "$dest_dir"
for dir in bin configs containers data workflows; do
  src="$src_dir/$dir/"
  dest="$dest_dir/$dir"

  if [ ! -d "$src" ]; then
    2>&1 echo "Folder '$src' does not exist!"
    exit 1
  fi

  ln -sf "$src" "$dest"
done
