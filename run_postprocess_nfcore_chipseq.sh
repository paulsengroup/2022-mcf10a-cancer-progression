#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x


step='postprocess_nfcore_chipseq'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

# shellcheck disable=SC2064
trap "rm -f '$wd/configs/postprocess_nfcore_chipseq.config'" EXIT

for config in configs/nfcore/postprocess_nfcore_chipseq_*.config; do
  ln -sf "$(readlink -f "$config")" "$wd/configs/postprocess_nfcore_chipseq.config"
  1>&2 echo "Processing $(basename "$config" .config)..."
  ./run_workflow.sh "$wd" "$step"
done
