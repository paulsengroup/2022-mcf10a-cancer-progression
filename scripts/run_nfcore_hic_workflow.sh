#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

argc="$#"

if [ $argc -lt 4 ]; then
  echo "Usage: $0 input_dir output_dir sample_name input_pattern"
  echo "Example: $0 2022-david-hic/ 2022-david-hic/data/scratch/nfcore-hic 10A HiC_001_10A{_R1,_R2}.fastq.gz --max_cpus=50 --max_memory=500.GB --max_time=300.h"
  exit 1
fi

indir="$1"
outdir="$2"
sample_name="$3"
input_name="$4"
args=("${@:5}")
data_dir="$indir/data"

wd="$NXF_WORK/nfcore_hic_$sample_name"
launch_dir="$indir/.launch-dir/$sample_name"
sentinel="$launch_dir/done"

mkdir -p "$wd" "$outdir/$sample_name" "$launch_dir"

if [ -f "$sentinel" ]; then
  echo 2>&1 "\"$sample_name\" has already been processed. SKIPPING!"
  exit 0
fi

scratch_dir="$(mktemp -d "$launch_dir/scratch-XXXXXX")"

function cleanup {
  exit_code=$?
  if [ $exit_code -ne 0 ]; then
    rm -f "$sentinel"
  fi

  rm -rf "$scratch_dir"
}

trap 'cleanup' EXIT

echo 2>&1 "Extracting bowtie2 index..." &&
  zstd -dcf "$data_dir/output/preprocessing/bowtie2_idx/GRCh38_genome_assembly.tar.zst" |
  tar -xf - -C "$scratch_dir/" &&
  echo 2>&1 "DONE"

bw2_idx_files=("$scratch_dir"*.b2)

bw2_idx="${bw2_idx_files[0]}"

rm -rf "$launch_dir/hic"
git clone --depth=1 --branch 20220513 https://github.com/robomics/hic.git "$launch_dir/hic"
ln -sf "$(readlink -f ./configs)/" "$launch_dir/configs"
ln -sf "$(readlink -f ./containers)/" "$launch_dir/containers"
ln -sf "$(readlink -f ./data)/" "$launch_dir/data"

# I am not sure who's creating these symlinks
sleep 5 && (unlink ./configs/configs || true) && (unlink ./data/data || true) && (unlink ./containers/containers || true) &

# Note: It is important that files are named like *_R[12].fastq.gz.
# Furthermore, file names shall not contain single digit fields:
# For example, mysample_001_R1.fastq.gz is ok, but mysample_1_R1.fastq.gz is not ok!
cd "$launch_dir"
nextflow run ./hic \
  "${args[@]}" \
  --input "$input_name" \
  --outdir "$outdir/$sample_name" \
  --bwt2_index "$bw2_idx" \
  -e."TMPDIR=$TMPDIR" \
  -c "$indir/configs/base_saga.config" \
  -c "$indir/configs/nfcore_hic.config" \
  -ansi-log false \
  -resume |&
  tee "$outdir/$sample_name/nfcore_hic.log"

echo 2>&1 "Touching \"$sentinel\"..."
touch "$sentinel"

nextflow clean -f $(nextflow log -q | tail -n 1)
