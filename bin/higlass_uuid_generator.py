#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import hashlib
import os
import sys
import uuid

# Generate reproducible UUIDs using hashed file name(s) as seed(s)
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(f"Usage: {sys.argv[0]} file1.txt ...")
        sys.exit(1)

    for f in sys.argv[1:]:
        hasher = hashlib.md5()
        hasher.update(os.path.basename(f).encode("utf-8"))

        print(uuid.UUID(hasher.hexdigest()))
