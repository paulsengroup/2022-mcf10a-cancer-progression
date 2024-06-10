#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u


wd=".nextflow-robomics-call-tad-cliques-wd1"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-v0.4.0.tar.xz' \
   configs/call_tad_cliques.json

wd=".nextflow-robomics-call-tad-cliques-wd2"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-v0.4.0.tar.xz' \
   configs/call_tad_cliques_vs_control.json


step='postprocess_call_tad_cliques'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running postprocessing step for robomics/call_tad_cliques...'
./run_workflow.sh "$wd" "$step"
