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
  'workflows/nfcore-rnaseq-v3.10.1.tar.xz' \
   configs/nfcore_rnaseq*.config
