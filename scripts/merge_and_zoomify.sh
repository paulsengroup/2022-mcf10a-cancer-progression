#!/usr/bin/env bash

set -e
set -u

argc="$#"

if [ $argc -lt 5 ]; then
  echo "Usage: $0 singularity.img numcpus resolution(s) output_path input_matrice(s)"
  echo "Example: $0 containers/cache/ghcr.io-paulsengroup-2022-david-hic-cooltools-0.5.1.img 8 1000N data/output/by_sample/HiC_10A/HiC_10A.mcool data/output/by_sample/HiC_???_10A/HiC_???_10A_raw.cool"
  exit 1
fi

container="$1"
ncpus="$2"
res="$3"
output_path="$4"
input_matrices=("${@:5}")

# Setup scratch space
SCRATCH="${SCRATCH-$TMPDIR}"
TMPDIR="${LOCALSCRATCH-$SCRATCH}"
wd="$(mktemp -d "$TMPDIR/wd-XXXXXXXXXX")"
# shellcheck disable=SC2064
trap "rm -rf '$wd'" EXIT
mkdir -p "$wd/input" "$wd/output"

# Merge coolers when appropriate, otherwise just copy the input cooler to the working dir
if [[ "${#input_matrices[@]}" -gt 1 ]]; then
  cp "${input_matrices[@]}" "$wd/input/"
  singularity exec -B "$wd:/data" "$container" \
    bash -c 'cooler merge /data/output/merged.cool /data/input/*.cool'
else
  cp "${input_matrices[0]}" "$wd/output/merged.cool"
fi


# Zoomify and balance
singularity exec -B "$wd:/data" "$container"  \
  cooler zoomify --nproc "$ncpus"             \
                 --resolutions "$res"         \
                 --balance                    \
                 --balance-args="-p $ncpus"   \
                 --out /data/output/out.mcool \
                 /data/output/merged.cool

# Copy output .mcool
mkdir -p "$(dirname "$output_path")"
cp "$wd/output/out.mcool"  "$output_path"
