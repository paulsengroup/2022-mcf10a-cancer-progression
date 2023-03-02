#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


wd=".nextflow-robomics-call-tad-cliques-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-v0.3.0.tar.xz' \
   configs/call_tad_cliques.config \
   configs/call_tad_cliques_vs_control.config
