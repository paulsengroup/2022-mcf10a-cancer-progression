#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -x
set -u

git_root="$(git rev-parse --show-toplevel)"

wd="$PWD"

trap 'cd "$wd"' EXIT

destdir="$git_root/data/output/nfcore_hic/mcools"
mkdir -p "$destdir"

cd "$destdir" || exit
for src in ../HiC*/*.mcool; do
  dest="$(basename "$src")"
  dest="GRCh38_${dest#HiC_}"
  ln -s "$src" "$dest"
done
