#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -x
set -o pipefail

trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

nextflow -version

mkdir -p logs/

function launch_runner_script {
  script_name="$1"
  log_name="${script_name#run_}"
  log_file="logs/${log_name%.sh}.log"

  2>&1 echo "Running \"$script_name\"..." && "./$script_name.sh" |& tee "$log_file" && 2>&1 echo "\"$script_name\" completed successfully!" || 2>&1 echo "\"$script_name\" failed!"
}

function prepare_raw_data {
  launch_runner_script run_fetch_data
  launch_runner_script run_preprocessing
}

function process_raw_data {
  launch_runner_script run_nfcore_rnaseq &
  launch_runner_script run_nfcore_hic && launch_runner_script run_postprocess_nfcore_hic &
  # TODO run chip workflow

  wait
}

function run_compartment_analysis {
  launch_runner_script run_compartment_analysis
}

function run_tad_analysis {
  launch_runner_script run_tad_analysis
  launch_runner_script run_call_tad_cliques_workflow
}

function run_de_analysis {
  launch_runner_script run_diff_expression_analysis
}

function run_comparative_analysis {
  launch_runner_script run_compartment_analysis
}

function run_chrom3d {
  launch_runner_script run_chrom3d
}

prepare_raw_data
# process_raw_data
run_tad_analysis
run_compartment_analysis
run_de_analysis
run_chrom3d
run_comparative_analysis
