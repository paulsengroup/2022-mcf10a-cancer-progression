#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo 1>&2 'Running nfcore/chipseq...'

wd=".nextflow-nfcore-chipseq-wd"
mkdir -p "$wd"

for dir in configs containers data scripts workflows; do
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
nextflow run nf-core/chipseq -r 2.0.0 \
  "${args[@]}" \
  -c configs/nfcore_chipseq_h3k27ac.config \
  -profile singularity \
  -resume
)

./remove_symlink_loops.sh
(cd "$wd" &&
nextflow run nf-core/chipseq -r 2.0.0 \
  "${args[@]}" \
  -c configs/nfcore_chipseq_h3kxxme.config \
  -profile singularity \
  -resume
)

