#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo 1>&2 'Running nfcore/rnaseq...'

wd=".nextflow-nfcore-rnaseq-wd"
mkdir -p "$wd"

(cd "$wd" && ln -sf ../configs/ configs)
(cd "$wd" && ln -sf ../containers/ containers)
(cd "$wd" && ln -sf ../data/ data)
(cd "$wd" && ln -sf ../scripts/ scripts)
(cd "$wd" && ln -sf ../workflows/ workflows)

if [[ $HOSTNAME == *.saga* ]]; then
    args=("${@:2}"
        --max_memory=400.G
        --max_cpus=52
        --max_time=336.h
        --project="${SLURM_PROJECT_ID-changeme}")
else
    args=()
fi

(cd "$wd" &&
nextflow run nf-core/rnaseq -r 3.9 \
  "${args[@]}" \
  --input data/nfcore_rnaseq_samplesheet.csv \
  --outdir data/output/nfcore_rnaseq/ \
  --aligner star_salmon \
  --pseudo_aligner salmon \
  --genome GRCh38 \
  -profile singularity \
  -resume
)
