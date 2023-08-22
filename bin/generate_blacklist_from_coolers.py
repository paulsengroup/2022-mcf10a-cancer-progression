#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import bioframe as bf
import cooler
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser("Output the list of blacklisted bins for one or multiple coolers.")

    cli.add_argument(
        "coolers",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more a .cool file (URI syntax supported).",
    )
    cli.add_argument(
        "--norm-name",
        type=str,
        default="weight",
        help="Name of the normalization to be applied.",
    )
    cli.add_argument(
        "--external-list",
        type=pathlib.Path,
        nargs="+",
        help="One or more BED files with intervals to include in the blacklist.",
    )

    return cli


def read_weights(coolers, norm_name) -> pd.DataFrame:
    series = []
    for clr in coolers:
        clr = cooler.Cooler(str(clr))
        bins = clr.bins()[:].set_index(["chrom", "start", "end"])
        series.append(bins[norm_name])

    return pd.concat(series, join="outer", axis="columns")


def main():
    args = vars(make_cli().parse_args())

    df = read_weights(args["coolers"], args["norm_name"])

    # Keep bins with NaN as weight in at least one cooler
    df = df[df.isna().any(axis="columns")].index.to_frame()

    # Read external blacklist
    if args["external_list"] is not None:
        dfs = [df]
        for path in args["external_list"]:
            dfs.append(bf.read_table(path, schema="bed3"))

        df = pd.concat(dfs)

    # Merge overlapping intervals
    df = bf.merge(df).drop(columns="n_intervals")

    df.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
