#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import re
import sys


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "chrom-bed",
        type=str,
        help="BED4+ file with the current and new chromosome names stored in columns 4 and 1 respectively",
    )

    return cli


def import_chrom_name_mappings(path_to_bed):
    mappings = {}
    with open(path_to_bed, "r") as f:
        for i, line in enumerate(f):
            toks = line.split("\t", 5)
            if len(toks) < 4:
                raise RuntimeError(
                    f'Line {i} from file "{path_to_bed}" does not seem to be in BED4+ format: '
                    f"expected 4 or more fields, found {len(toks)}"
                )

            id1 = toks[3].strip()
            id2 = toks[0].strip()
            if id1 in mappings:
                raise RuntimeError(
                    f'Found a duplicate entry for "{id1}" at line {i} of file "{path_to_bed}"'
                )

            mappings[id1] = id2

    if len(mappings) == 0:
        raise RuntimeError(f'Unable to import any chromosome from file "{path_to_bed}"')

    return mappings


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    pattern = re.compile(r"^>(\S+).*$")

    chrom_ids = import_chrom_name_mappings(args["chrom-bed"])

    for line in sys.stdin:
        match = pattern.match(line)
        if match:
            old_id = match.group(1)
            if old_id in chrom_ids:
                new_id = chrom_ids[old_id]
                line = line.replace(old_id, new_id, 1)
        print(line, end="")
