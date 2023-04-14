#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

wd=".nextflow-nfcore-hic-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

echo 1>&2 'Running nfcore/hic...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-nfcore-hic-v2.0.0-patched.tar.xz' \
  configs/nfcore_hic*.config
