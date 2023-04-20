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

cd "$destdir" || exit 1
for src in ../HiC*/*.mcool ../*merged/*.mcool; do
  dest="$(basename "$src")"
  dest="hg38_${dest#HiC_}"
  ln -sf "$src" "$dest"
done

destdir="$git_root/data/output/nfcore_hic/hic"
mkdir -p "$destdir"

cd "$destdir" || exit 1
for src in ../HiC*/*.hic ../*merged/*.hic; do
  dest="$(basename "$src")"
  dest="hg38_${dest#HiC_}"
  ln -sf "$src" "$dest"
done
