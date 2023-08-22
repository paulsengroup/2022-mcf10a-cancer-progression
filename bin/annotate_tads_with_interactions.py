#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings
from typing import Union

import bioframe as bf
import cooler
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "cooler",
        type=pathlib.Path,
        help="Path to a .cool file (URI syntax supported).",
    )

    cli.add_argument(
        "tads",
        type=pathlib.Path,
        help="Path to TADs.",
    )

    return cli


def import_bed(path_to_bed: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path_to_bed, schema="bed3")


def get_interactions(coords, clr: cooler.Cooler) -> Union[int, float]:
    chrom, start, end = coords
    return clr.matrix(balance=False, sparse=True).fetch(f"{chrom}:{start}-{end}").sum()


def main():
    args = vars(make_cli().parse_args())

    tads = import_bed(args["tads"])
    clr = cooler.Cooler(str(args["cooler"]))

    tads["value"] = tads.apply(get_interactions, args=(clr,), axis="columns")

    tads.to_csv(sys.stdout, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
