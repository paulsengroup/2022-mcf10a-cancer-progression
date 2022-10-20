#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x

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

    ln -sf ../configs/ configs
    ln -sf ../containers/ containers
    ln -sf ../data/ data
    ln -sf ../scripts/ scripts
    ln -sf ../workflows/ workflows

    nextflow run \
        "${args[@]}" \
        -c "configs/$name.config" \
        -c "$base_config" \
        -process.cache=deep \
        "workflows/$name.nf" \
        -ansi-log \
        -resume
}

step='compartment_analysis'
echo 1>&2 "Running step $step..."

wd=".nextflow-$step-wd"
mkdir -p "$wd"
(cd "$wd" && run_workflow "$step")
echo 1>&2 "Workflow \"$step\" successfully completed!"
