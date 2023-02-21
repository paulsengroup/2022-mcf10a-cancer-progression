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
    args=("${@:2}"
        --max_memory=400.G
        --max_cpus=52
        --max_time=336.h
        --project="${SLURM_PROJECT_ID-changeme}")
else
    args=()
fi

./remove_symlink_loops.sh
(cd "$wd" &&
nextflow run https://github.com/robomics/call_tad_cliques \
  -r v0.1.0 \
  "${args[@]}" \
  -c configs/call_tad_cliques_chrom3d.config \
  -resume
)

exit 0  # TODO: removeme!

function run_workflow() {
    name="$1"
    if [[ $HOSTNAME == *.saga* ]]; then
        base_config='configs/base_saga.config'
        args=("${@:2}"
            --max_memory=400.G
            --max_cpus=52
            --max_time=336.h
            --project="${SLURM_PROJECT_ID-changeme}")
    else
        base_config='configs/base_hovig.config'
        args=("${@:2}")
    fi

    for dir in bin configs containers data workflows; do
        ln -sf "../$dir/" "$dir"
    done

    ../remove_symlink_loops.sh
    nextflow run \
        "${args[@]}" \
        -c "configs/$name.config" \
        -c "$base_config" \
        -process.cache=deep \
        "workflows/$name.nf" \
        -ansi-log \
        -resume
}

echo 1>&2 "Running step $step..."
(cd "$wd" && run_workflow "$step")
echo 1>&2 "Workflow \"$step\" successfully completed!"
