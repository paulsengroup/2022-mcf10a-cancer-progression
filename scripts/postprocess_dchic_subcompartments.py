#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedgraph",
        type=pathlib.Path,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
    )

    return cli


def import_data(path_to_df: pathlib.Path) -> pd.DataFrame:
    if str(path_to_df) == "-":
        path_to_df = sys.stdin

    return pd.read_table(path_to_df).rename(columns={"chr": "chrom"})


def compute_state_mode(df: pd.DataFrame) -> pd.Series:
    mode = df.filter(regex=".state$").mode(axis="columns")
    num_conditions = len(mode.columns)

    # df.mode() returns a df with the same number of columns as the input df.
    # Columns 2+ are used to report ties.
    # When there are no ties (or at the very least, fewer ties than there are columns),
    # values are set to nan
    mask = mode.isna().sum(axis="columns") != (num_conditions - 1)

    mode.loc[mask, 0] = "None"
    return mode[0]


def main():
    args = vars(make_cli().parse_args())

    df = import_data(args["bedgraph"])

    df["state.mode"] = compute_state_mode(df)
    df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
