#!/usr/bin/env python3

# Copyright (C) 2023 Jonas Paulsen <jonaspau@uio.no>
#
# SPDX-License-Identifier: MIT

import re
import sys

gtrack_fname = str(sys.argv[1])

with open(gtrack_fname) as gf:
    gtrack_file = gf.readlines()

for line in gtrack_file:
    line = line.rstrip()
    if line.startswith("#"):
        print(line)
    else:
        line1 = re.sub(r"(chr\w+)", r"\1_A", line)
        line2 = re.sub(r"(chr\w+)", r"\1_B", line)
        print(line1)
        print(line2)
