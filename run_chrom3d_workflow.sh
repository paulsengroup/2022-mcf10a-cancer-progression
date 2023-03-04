#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


step='chrom3d'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-v0.3.0.tar.xz' \
   configs/call_tad_cliques_chrom3d.config


./run_workflow.sh "$wd" "$step"
