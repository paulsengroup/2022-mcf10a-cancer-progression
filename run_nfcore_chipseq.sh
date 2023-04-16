#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


wd=".nextflow-nfcore-chipseq-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running nfcore/chipseq...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/nfcore-chipseq-v2.0.0.tar.xz' \
  configs/nfcore/nfcore_chipseq*.json
