#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"
sudo -E higlass-manage stop "$name"
