#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


wd=".nextflow-robomics-chrom3d-nf-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

echo 1>&2 'Running robomics/chrom3d-nf...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-chrom3d-nf-v0.0.2.tar.xz' \
  configs/chrom3d.json
