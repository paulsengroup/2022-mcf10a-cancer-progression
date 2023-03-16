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

bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP1_WT MCF10A_WT_REP1 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP2_WT MCF10A_WT_REP2 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP3_WT MCF10A_WT_REP3 --replace

bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP1_T1 MCF10A_T1_REP1 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP2_T1 MCF10A_T1_REP2 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP3_T1 MCF10A_T1_REP3 --replace

bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP1_C1 MCF10A_C1_REP1 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP2_C1 MCF10A_C1_REP2 --replace
bin/fix_nfcore_rnaseq_output.sh "$wd/data/output/nfcore_rnaseq" MCF10A_REP3_C1 MCF10A_C1_REP3 --replace
