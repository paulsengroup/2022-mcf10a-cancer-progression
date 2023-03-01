#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


step='chrom3d'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

echo 1>&2 'Running robomics/call_tad_cliques...'

for dir in bin configs containers data workflows; do
  (cd "$wd" && ln -sf "../$dir/" "$dir")
done

if [[ $HOSTNAME == *.saga* ]]; then
    base_config='configs/base_saga.config'
    args=("${@:2}"
        --max_memory=400.G
        --max_cpus=52
        --max_time=336.h
        --project="${SLURM_PROJECT_ID-changeme}")
else
    base_config='configs/base_hovig.config'
    args=()
fi

./remove_symlink_loops.sh
(cd "$wd" &&
nextflow run https://github.com/robomics/call_tad_cliques \
  -r v0.3.0 \
  "${args[@]+"${args[@]}"}" \
  -c "$base_config" \
  -c configs/call_tad_cliques_chrom3d.config \
  -resume
)

exit 0  # TODO: removeme!

./setup_workflow_workdir.sh "$PWD" "$wd"
./run_workflow.sh "$wd" "$step"
