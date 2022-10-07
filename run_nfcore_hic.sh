#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

function run_workflow() {
  set -e
  set -u
  set -o pipefail

  fq="$1"

  bname="${fq%_R1.fastq.gz}"
  sample_name="$(basename "$bname")"
  outdir="$(readlink -f ./data/scratch/nfcore_hic)"
  input_name="${bname}{_R1,_R2}.fastq.gz"

  printf 2>&1 "Processing %s...\n" "$sample_name"

  if [[ $HOSTNAME == *.saga* ]]; then
    sbatch -A "${SLURM_PROJECT_ID-changeme}" \
      -c 2 --mem 3G -t 4-00:00:00 \
      --job-name="nfcore_hic_$sample_name" \
      --wait \
      "scripts/run_nfcore_hic_workflow.sh" \
      "$(readlink -f .)" \
      "$outdir" \
      "$sample_name" \
      "$input_name" \
      --process.cache=deep \
      --max_cpus 52 \
      --max_memory '400.GB' \
      --max_time 336.h &>/dev/null
  else
    scripts/run_nfcore_hic_workflow.sh \
      "$(readlink -f .)" \
      "$outdir" \
      "$sample_name" \
      "$input_name" \
      --process.cache=deep \
      --max_cpus 128 \
      --max_memory '500.GB' \
      --max_time 336.h &>/dev/null
  fi

  printf 2>&1 "DONE processing %s!\n" "$sample_name"
}

export -f run_workflow

if [[ $HOSTNAME == *.saga* ]]; then
  parallel_runs=999
else
  parallel_runs=2
fi

mkdir -p ./data/scratch/nfcore_hic
raw_reads_dir="$(readlink -f ./data/input/raw_data/hic)"
# Note: It is important that files are named like *_R[12].fastq.gz.
# Furthermore, file names shall not contain single digit fields:
# For example, mysample_001_R1.fastq.gz is ok, but mysample_1_R1.fastq.gz is not ok!
printf "%s\n" "$raw_reads_dir/"*R1.fastq.gz |
  xargs -n 1 -P "$parallel_runs" -I {} bash -c 'run_workflow "$@" || exit 255' _ {}
