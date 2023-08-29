#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings
from typing import List

import bioframe as bf
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "insulation",
        nargs="+",
        type=pathlib.Path,
        help="Path to two or more bedgraph with the scores to be processed.",
    )

    cli.add_argument(
        "--blacklist",
        type=pathlib.Path,
        help="Path to a BED file with the list of regions to be masked",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma separated list of labels to use instead of file names.",
    )

    return cli


def import_bedgraph(path_to_bedgraph: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path_to_bedgraph, schema="bedGraph").set_index(["chrom", "start", "end"])


def import_insulation_scores(paths: List[pathlib.Path], labels: List[str]) -> pd.DataFrame:
    assert len(paths) == len(labels)
    dfs = []
    for label, path in zip(labels, paths):
        dfs.append(import_bedgraph(path).rename(columns={"value": label}))

    return pd.concat(dfs, axis="columns", join="outer")


def main():
    args = vars(make_cli().parse_args())

    labels = args.get("labels")
    if labels is not None:
        labels = labels.split(",")
    paths_to_scores = args["insulation"]

    if labels is None:
        labels = [p.name for p in paths_to_scores]
    else:
        if len(labels) != len(paths_to_scores):
            raise RuntimeError(
                f"Mismatch in the number of files and labels: expected {len(paths_to_scores)} labels, found {len(labels)}"
            )

    scores = import_insulation_scores(paths_to_scores, labels)

    if args["blacklist"] is not None:
        blist = bf.read_table(args["blacklist"], schema="bed3")
        scores = bf.setdiff(scores.reset_index(), blist).set_index(["chrom", "start", "end"])

    scores["std"] = scores.std(axis="columns")
    scores.to_csv(sys.stdout, sep="\t", index=True, header=True, na_rep="NaN")


if __name__ == "__main__":
    main()
