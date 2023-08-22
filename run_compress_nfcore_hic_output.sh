#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

wd=".nextflow-nfcore-hic-output-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

echo 1>&2 'Running robomics/compress-nfcore-hic-output...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-compress-nfcore-hic-output-v0.0.1.tar.xz' \
  configs/compress_nfcore_hic_output.json
