#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo 1>&2 'Running nfcore/nascent...'

wd=".nextflow-nfcore-nascent-wd"
mkdir -p "$wd"

for dir in bin configs containers data workflows; do
    (cd "$wd" && ln -sf "../$dir/" "$dir")
done

if [[ $HOSTNAME == *.saga* ]]; then
    args=("${@:2}"
        --max_memory=400.GB
        --max_cpus=52
        --max_time=336.h
        --project="${SLURM_PROJECT_ID-changeme}")
else
    args=()
fi

./remove_symlink_loops.sh
(cd "$wd" &&
nextflow run nf-core/nascent -r 2.1.1 \
  "${args[@]}" \
  -c configs/nfcore_nascent.config \
  -profile singularity \
  -resume
)