#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

function run_workflow() {
    name="$1"
    base_config='configs/base_saga.config'
    args=("${@:2}"
          --max_memory=400.G
          --max_cpus=52
          --max_time=336.h
          --project="${SLURM_PROJECT_ID-changeme}")

    nextflow run \
        "${args[@]}" \
        -c "configs/$name.config" \
        -c "$base_config" \
        "workflows/$name.nf" \
        -resume
}


preproc_steps=(fetch_data
               preprocess_data
)

for step in "${preproc_steps[@]}"; do
  run_workflow "$step"
done


mkdir -p ./data/output/nfcore_hic
raw_reads_dir="$(readlink -f ./data/input/raw_data/hic)"
# Note: It is important that files are named like *_R[12].fastq.gz.
# Furthermore, file names shall not contain single digit fields:
# For example, mysample_001_R1.fastq.gz is ok, but mysample_1_R1.fastq.gz is not ok!
for fq in "$raw_reads_dir/"*R1.fastq.gz; do
  bname="${fq%_R1.fastq.gz}"
  sample_name="$(basename "$bname")"
  outdir="$(readlink -f ./data/output/nfcore_hic)"
  input_name="${bname}_R{1,2}.fastq.gz"

  sbatch -A "${SLURM_PROJECT_ID-changeme}" \
    -c 2 --mem 3G -t 4-00:00:00            \
    --job-name="nfcore_hic_$sample_name"   \
    --wait                                 \
    "scripts/run_nfcore_hic.sh"            \
    "$(readlink -f .)"                     \
    "$outdir"                              \
    "$sample_name"                         \
    "$input_name" &
done
wait

scripts/symlink_nfcore_hic_output.sh "$(readlink -f ./data/output/nfcore_hic)" \
                                     "data/output/by_sample"
