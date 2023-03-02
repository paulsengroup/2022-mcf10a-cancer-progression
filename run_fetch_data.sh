#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x


step='fetch_data'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"
./run_workflow.sh "$wd" "$step"
