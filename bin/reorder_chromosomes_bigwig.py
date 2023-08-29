#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings
from typing import List, Tuple, Union

import bioframe as bf
import cooler
import pandas as pd
import pyBigWig


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bigwig",
        type=pathlib.Path,
        help="Path to the .bigwig file to be re-ordered.",
    )

    cli.add_argument(
        "chrom-sizes",
        type=pathlib.Path,
        help="Path to a .chrom.sizes file.",
    )

    cli.add_argument(
        "-o",
        "--output",
        required=True,
        type=pathlib.Path,
        help="Output path (including the extension)",
    )

    return cli


def import_chromsizes(path_to_chrom_sizes: pathlib.Path) -> List[Tuple[str, int]]:
    df = pd.read_table(path_to_chrom_sizes, names=["chrom", "length"]).set_index("chrom")
    chroms = []
    for chrom, size in df.itertuples():
        chroms.append((chrom, size))

    return chroms


def main():
    args = vars(make_cli().parse_args())

    chroms = import_chromsizes(args["chrom-sizes"])

    with pyBigWig.open(str(args["bigwig"])) as bwi, pyBigWig.open(str(args["output"]), "w") as bwo:
        bwo.addHeader(chroms)
        for chrom, length in chroms:
            if chrom not in bwi.chroms():
                continue

            starts = []
            ends = []
            values = []
            for start, end, value in bwi.intervals(chrom, 0, length):
                starts.append(start)
                ends.append(end)
                values.append(value)
            chroms_ = len(starts) * [chrom]

            bwo.addEntries(chroms_, starts, ends=ends, values=values)


if __name__ == "__main__":
    main()
