#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

argc="$#"

if [ $argc -ne 4 ]; then
  echo "Usage: $0 input_dir output_dir sample_name input_pattern"
  echo "Example: $0 2022-david-hic/ 2022-david-hic/data/output/nfcore-hic 10A HiC_001_10A_R{1,2}.fq.gz"
  exit 1
fi

indir="$1"
outdir="$2"
sample_name="$3"
input_name="$4"
data_dir="$indir/data"
bw2_idx="$data_dir/output/preprocessing/bowtie2_idx/GRCh38_genome_assembly"

wd="$NXF_WORK/nfcore_hic_$sample_name"
launch_dir="$indir/.launch-dir/$sample_name"
mkdir -p "$wd" "$outdir/$sample_name" "$launch_dir"

rsync -aPq "$HOME/.nextflow/assets/robomics/hic" "$launch_dir/" --delete
(cd "$launch_dir/hic" && git checkout 20220513)
ln -sf "$(readlink -f ./containers)/" "$launch_dir/containers"
ln -sf "$(readlink -f ./data)/" "$launch_dir/data"

# Note: It is important that files are named like *_R[12].fastq.gz.
# Furthermore, file names shall not contain single digit fields:
# For example, mysample_001_R1.fastq.gz is ok, but mysample_1_R1.fastq.gz is not ok!
cd "$launch_dir"
nextflow run ./hic                             \
         --max_cpus 52                         \
         --max_memory '400.GB'                 \
         --max_time 336.h                      \
         --input "$input_name"                 \
         --outdir "$outdir/$sample_name"       \
         --bwt2_index "$bw2_idx"               \
         -e."TMPDIR=$TMPDIR"                   \
         -c "$indir/configs/base_saga.config"  \
         -c "$indir/configs/nfcore_hic.config" \
         -ansi-log false                       \
         -resume |&
         tee "$outdir/$sample_name/nfcore_hic.log"
# Remove write permissions
chmod -R a-w "$outdir/$sample_name"

# rm -r "$launch_dir"
