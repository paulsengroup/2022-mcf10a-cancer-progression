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

# shellcheck disable=SC2207
IFS=$'\n' \
sample_names=($(echo ./data/output/nfcore_hic/* |
                tr ' ' '\n' |
                sed -E 's/.*HiC_[[:digit:]]{3}_(.*)$/\1/g' |
                sort -u)
             )

function run_merge_and_zoomify {
  cooler_image='containers/cache/ghcr.io-paulsengroup-2022-david-hic-cooltools-0.5.1.img'
  outname="$1"
  cooler=("${@:2}")

  job_name="zoomify_and_merge_$(basename "$outname" .mcool)"

  sbatch -A "${SLURM_PROJECT_ID-changeme}" \
    -c 8 --mem 20G -t 04:00:00             \
    --gres=localscratch:15G                \
    --job-name="$job_name"                 \
    --wait                                 \
    "scripts/merge_and_zoomify.sh"         \
    "$cooler_image"                        \
    8 1000N "$outname" "${cooler[@]}"
}

for sample_name in "${sample_names[@]}"; do
  coolers=(./data/output/by_sample/HiC_???_"$sample_name"/*_raw.cool)
  for cooler in "${coolers[@]}"; do

    outname="${cooler%_raw.cool}.mcool"
    if [ ! -f "$outname" ]; then
      cooler=("$cooler")
      run_merge_and_zoomify "$outname" "${cooler[@]}" &
    fi
  done

  outname="data/output/by_sample/HiC_$sample_name/HiC_$sample_name.mcool"
  if [ ! -f "$outname" ]; then
    run_merge_and_zoomify "$outname" "${coolers[@]}" &
  fi
done
wait
