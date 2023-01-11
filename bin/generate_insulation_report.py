#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import warnings
from math import inf
from typing import List

import bioframe as bf
import matplotlib.pyplot as plt
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedgraph",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more bedgraph file(s) with the insulation scores to be processed.",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma separated list of sample labels to use instead of those inferred from file names.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
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


def detect_path_collisions(*args: pathlib.Path) -> None:
    for path in args:
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite file {path}. Pass --force to overwrite existing file(s).")


def main():
    args = vars(make_cli().parse_args())

    output_prefix = args["output_prefix"]
    if not args["force"]:
        detect_path_collisions(output_prefix.with_suffix(".png"), output_prefix.with_suffix(".svg"))

    labels = args.get("labels")
    paths_to_scores = args["bedgraph"]
    if labels is None:
        labels = [p.name for p in paths_to_scores]
    else:
        labels = labels.split(",")
        if len(labels) != len(paths_to_scores):
            raise RuntimeError(
                f"Mismatch in the number of files and labels: expected {len(paths_to_scores)} labels, found {len(labels)}"
            )

    scores = import_insulation_scores(paths_to_scores, labels)

    grid_size = len(scores.columns)

    fig, axs = plt.subplots(grid_size, grid_size, figsize=(grid_size * 6.4, grid_size * 6.4))

    lb = inf
    ub = -inf
    for i in range(grid_size):
        for j in range(i + 1, grid_size):
            cond1, cond2 = scores.columns[i], scores.columns[j]
            df = scores[[cond1, cond2]].dropna()
            axs[i][j].scatter(df[cond1], df[cond2])

            lb = min(lb, df[cond1].min(), df[cond2].min())
            ub = max(ub, df[cond1].max(), df[cond2].max())

    ub += abs(ub) * 0.05
    lb -= abs(lb) * 0.05
    for i in range(grid_size):
        for j in range(grid_size):
            ax = axs[i][j]

            if i == 0:
                ax.set(title=scores.columns[j])

            if j == i + 1:
                ax.set(ylabel=scores.columns[i])

            if i >= j:
                fig.delaxes(ax)
            else:
                ax.set(xlim=(lb, ub), ylim=(lb, ub))

    plt.tight_layout()

    fig.savefig(output_prefix.with_suffix(".png"), dpi=600)
    fig.savefig(output_prefix.with_suffix(".svg"))


if __name__ == "__main__":
    main()
