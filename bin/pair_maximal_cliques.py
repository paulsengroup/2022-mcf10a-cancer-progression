#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import pathlib
import sys
from typing import List

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        nargs="+",
        help="Path to the TAD clique TSV (cliques should be called using same TAD annotation).",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma-separated list of labels used to refer to files. Should be in the same order as the provided files.",
    )

    return cli


def import_data(path_to_dfs: List, labels: List) -> List[pd.DataFrame]:
    dfs = []

    for path, label in zip(path_to_dfs, labels):
        df = pd.read_table(path).set_index("name").rename(columns={"size": label})
        df["tad_ids"] = df["tad_ids"].str.split(",")
        dfs.append(df)

    return dfs


def compute_max_tad_clique_size(cliques: pd.DataFrame, label: str) -> pd.DataFrame:
    clique_sizes = {}

    for _, (tads, size) in cliques[["tad_ids", label]].iterrows():
        for tad in tads:
            if tad in clique_sizes:
                clique_sizes[tad] = max(clique_sizes[tad], size)
            else:
                clique_sizes[tad] = size

    df = pd.DataFrame({"tad": clique_sizes.keys(), label: clique_sizes.values()})
    df["tad"] = df["tad"].astype(int)

    return df.set_index("tad").sort_index()


def main():
    args = vars(make_cli().parse_args())
    labels = args["labels"].split(",")
    dfs = import_data(args["tsv"], labels)

    dfs = [compute_max_tad_clique_size(df, label) for label, df in zip(labels, dfs)]

    df = functools.reduce(
        lambda left, right: pd.merge(left, right, how="left", left_index=True, right_index=True),
        dfs,
    )

    df.to_csv(sys.stdout, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()
