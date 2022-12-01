#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo 1>&2 'Running robomics/call_tad_cliques...'

wd=".nextflow-robomics-call-tad-cliques-wd"
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
nextflow run https://github.com/robomics/call_tad_cliques \
  -r v0.0.4 \
  "${args[@]}" \
  -c configs/call_tad_cliques.config \
  -resume
)

(cd "$wd" &&
nextflow run https://github.com/robomics/call_tad_cliques \
  -r v0.0.4 \
  "${args[@]}" \
  -c configs/call_tad_cliques_vs_control.config \
  -resume
)
