#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


wd=".nextflow-nfcore-nascent-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running nfcore/nascent...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/nfcore-nascent-v2.1.1.tar.xz' \
   configs/nfcore_nascent*.config
