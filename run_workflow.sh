#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

# IMPORTANT! This script should be run from the repository root!

if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 working-dir workflow-name"
  2>&1 echo "Example: $0 .nextflow-fetch_data-wd fetch_data"
  exit 1
fi

if ! command -v nextflow &> /dev/null; then
  2>&1 echo "Unable to find nextflow in your PATH"
  exit 1
fi

wd="$1"
name="$2"
workflow="workflows/$name.nf"
config="configs/$name.config"

if [ ! -d "$wd" ]; then
  2>&1 echo "Working directory '$wd' does not exist! Please create it before running $(basename "$0")"
  exit 1
fi

cwd="$PWD"
# shellcheck disable=SC2064
trap "cd '$cwd'" EXIT
cd "$wd"

if [ ! -f "$workflow" ]; then
  2>&1 echo "Unable to find workflow '$name'!"
  2>&1 echo "File '$workflow' does not exist"
  exit 1
fi

if [ ! -f "$config" ]; then
  2>&1 echo "Unable to find config for workflow '$name'!"
  2>&1 echo "File '$config' does not exist"
  exit 1
fi

args=("${@:3}")
if [[ $HOSTNAME == *.saga* ]]; then
  base_config="configs/base_saga.config"
  args+=(--max_memory=400.G
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
  2>&1 echo "Unable to find base config '$base_config': no such file"
  exit 1
fi

# This copy is important!
# Otherwise Nextflow won't add scripts located under bin/ to the runner PATH
cp "workflows/$name.nf" "$name.nf"

2>&1 echo "### Running step '$name'..."
cmd=(nextflow run
     "${args[@]+"${args[@]}"}"
     -c "configs/$name.config"
     -c "$base_config"
     -process.cache=deep
     "$name.nf"
     -ansi-log
     -resume)

"$cwd/remove_symlink_loops.sh"

2>&1 echo "${cmd[@]}"
"${cmd[@]}"

2>&1 echo "### Workflow \"$name\" successfully completed!"