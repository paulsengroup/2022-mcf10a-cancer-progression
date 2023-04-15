#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import re
import sys
from typing import List, Set, Tuple

import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("glob-pattern", type=str, help="Glob pattern matching files to be included in the samplesheet."),

    cli.add_argument(
        "--suffix", type=str, default=r"(\.part_\d+)?\.fastq.*$", help="Regex expression matching file name suffixes."
    ),
    cli.add_argument(
        "--mate-suffixes",
        type=str,
        default="R1,R2",
        help="Comma-separated list of suffixes used to discern between mate 1 and mate 2 FASTQ files.",
    )

    return cli


def generate_file_list(pattern: str) -> Set[pathlib.Path]:
    dirname = pathlib.Path(pattern).parent
    glob_pattern = pathlib.Path(pattern).name

    return set(dirname.glob(glob_pattern))


def group_mates(files: Set[pathlib.Path], mate_suffixes: str, file_name_suffix: str) -> List[Tuple[str, str, str]]:
    m1_suffix, _, m2_suffix = mate_suffixes.partition(",")
    suffix_pattern = re.compile(file_name_suffix)

    pairs = []
    for f in files:
        # Extract prefix and suffix (suffix looks something like .part_025.fastq.zst)
        suffix = suffix_pattern.search(str(f)).group(0)
        prefix = str(f).removesuffix(suffix)

        is_mate1 = prefix.endswith(m1_suffix)
        if is_mate1:
            prefix = prefix.removesuffix(m1_suffix)
            mate1 = f
            mate2 = prefix + m2_suffix + suffix

            assert pathlib.Path(mate2) in files, mate2
            pairs.append((prefix.strip("_"), mate1, mate2))

    return pairs


def main():
    args = vars(make_cli().parse_args())

    files = generate_file_list(args["glob-pattern"])

    paired_mates = group_mates(files, args["mate_suffixes"], args["suffix"])

    df = pd.DataFrame(paired_mates, columns=["sample", "fastq_1", "fastq_2"]).sort_values(
        by=["sample", "fastq_1", "fastq_2"]
    )

    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
