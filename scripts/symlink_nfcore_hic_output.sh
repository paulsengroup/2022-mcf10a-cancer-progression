#!/usr/bin/env bash

set -e
set -u
# set -x

argc="$#"

if [ $argc -ne 2 ]; then
  echo "Usage: $0 input_dir output_dir"
  echo "Example: $0 data/output/nfcore-hic data/output/samples"
  exit 1
fi

input_dir="$1"
output_dir="$2"

for dir in  "$input_dir/HiC_"*; do
  sample_name="$(basename "$dir")"

  mkdir -p "$output_dir/$sample_name"

  ln -sf "$dir/contact_maps/raw/cool/${sample_name}_R_1000.cool" \
         "$output_dir/$sample_name/${sample_name}_raw.cool"
done