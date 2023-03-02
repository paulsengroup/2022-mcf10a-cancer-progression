#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"

data_dir="$git_root/data/"

wd="$data_dir/output/tad_analysis/clodius"
rm -rf "$wd"
mkdir -p "$wd"

chrom_sizes="$data_dir/output/preprocessing/chrom_sizes/GRCh38.chrom.sizes"

clodius_wrapper="$(mktemp)"
trap 'rm -f "$clodius_wrapper"' EXIT

chmod 755 "$clodius_wrapper"

cat << 'EOF' >> "$clodius_wrapper"
#!/usr/bin/env bash

set -e
set -o pipefail
set -u

wd="$1"
chrom_sizes="$2"
bed="$3"
# bed is a path like tad_analysis/ICE/10000/*_domains.bed.gz

resolution="$(basename "$(dirname "$bed")")"
normalization="$(basename "$(dirname "$(dirname "$bed")")")"

out_name="$wd/$(basename "$bed" .bed.gz)_${normalization}_${resolution}.beddb"

echo "Processing \"$bed\"..."
clodius aggregate bedfile \
  --chromsizes-filename="$chrom_sizes" \
  -o "$out_name" "$bed" > /dev/null
EOF

function run_clodius {
   "$clodius_wrapper" "$wd" "$chrom_sizes" "$1"
}

export -f run_clodius

printf '%s\n' "$data_dir/output/tad_analysis/"*/*/*_domains.bed.gz |
xargs -L 1 -P "$(nproc)" bash -c "\"$clodius_wrapper\" \"$wd\" \"$chrom_sizes\" \"\$1\"" bash
