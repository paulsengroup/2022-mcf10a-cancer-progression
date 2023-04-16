#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

# IMPORTANT! This script should be run from the repository root!

if [ $# -lt 3 ]; then
  1>&2 echo "Usage: $0 working-dir workflow-archive param_files..."
  1>&2 echo "Example: $0 .nextflow-nfcore-chipseq-wd nfcore-chipseq-v2.0.0.tar.xz configs/nfcore/nfcore_chipseq_*.json"
  exit 1
fi

if ! command -v nextflow &> /dev/null; then
  1>&2 echo "Unable to find nextflow in your PATH"
  exit 1
fi

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}


wd="$1"
workflow_archive="$(readlink_py "$2")"
name="$(basename "$workflow_archive" .tar.xz)"
param_files=("${@:3}")

if [ ! -d "$wd" ]; then
  1>&2 echo "Working directory '$wd' does not exist! Please create it before running $(basename "$0")"
  exit 1
fi

if [ ! -f "$workflow_archive" ]; then
  1>&2 echo "Archive '$workflow_archive' does not exist!"
  exit 1
fi

for param_file in "${param_files[@]}"; do
  if [ ! -f "$param_file" ]; then
    1>&2 echo "Unable to find params file '$params_file'!"
    exit 1
  fi
done

cwd="$PWD"
# shellcheck disable=SC2064
trap "cd '$cwd'" EXIT
cd "$wd"

rm -rf "${name:?}/"
mkdir "$name/"
tar -xf "$workflow_archive" -C "$name/" --strip-components 1


args=()
if [[ $HOSTNAME == *.saga* ]]; then
  base_config="configs/base_saga.config"
  args+=(--max_memory=400.GB
         --max_cpus=52
         --max_time=336.h
         --project="${SLURM_PROJECT_ID-changeme}")
elif [[ $HOSTNAME == hovig*.uio.no ]]; then
  base_config="configs/base_hovig.config"
elif [[ "$OSTYPE" == darwin* ]]; then
  base_config="configs/base_macos.config"
else
  base_config="configs/base_linux.config"
fi

if [ ! -f "$base_config" ]; then
  1>&2 echo "Unable to find base config '$base_config': no such file"
  exit 1
fi

for params_file in "${param_files[@]}"; do
  1>&2 echo "### Running workflow '$name' -params-file '$(basename "$params_file")'..."
  cmd=(nextflow run
       "${args[@]+"${args[@]}"}"
       -params-file "$params_file"
       -c "$base_config"
       -process.cache=deep
       "$name"
       -ansi-log
       -resume)

  "$cwd/remove_symlink_loops.sh"

  1>&2 echo "${cmd[@]}"
  "${cmd[@]}"

  1>&2 echo "### Workflow \"$name\" successfully completed!"
done
