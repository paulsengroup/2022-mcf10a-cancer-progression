#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x


step='fish'
wd=".nextflow-$step-wd"
mkdir -p "$wd"

./setup_workflow_workdir.sh "$PWD" "$wd"
./run_workflow.sh "$wd" "$step"
