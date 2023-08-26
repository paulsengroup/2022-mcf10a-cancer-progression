#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

wd=".nextflow-nfcore-rnaseq-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running nfcore/rnaseq...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/nfcore-rnaseq-v3.11.2.tar.xz' \
  configs/nfcore/nfcore_rnaseq*.json

step='postprocess_nfcore_rnaseq'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

echo 1>&2 'Running postprocessing for nfcore/rnaseq...'
./run_workflow.sh "$wd" "$step"
