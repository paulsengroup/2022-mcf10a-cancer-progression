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
  'workflows/robomics-call_tad_cliques-0.5.0.tar.xz' \
   configs/call_tad_cliques_merged.json

wd=".nextflow-robomics-call-tad-cliques-wd2"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-0.5.0.tar.xz' \
   configs/call_tad_cliques_repl.json


wd=".nextflow-robomics-call-tad-cliques-wd3"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running robomics/call_tad_cliques...'
./run_external_workflow.sh \
  "$wd" \
  'workflows/robomics-call_tad_cliques-0.5.0.tar.xz' \
   configs/call_tad_cliques_merged_with_masking.json


step='postprocess_call_tad_cliques'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"

1>&2 echo 'Running postprocessing step for robomics/call_tad_cliques...'
./run_workflow.sh "$wd" "$step"
