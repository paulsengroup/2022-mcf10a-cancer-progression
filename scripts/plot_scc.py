#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to the TSV produced by run_hicrep.py.",
    )

    cli.add_argument(
        "--chrom-sizes",
        type=pathlib.Path,
        help="Path to a .chrom.sizes files. Used to compute the weighted average of correlation coefficients.",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        type=pathlib.Path,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def import_tsv(path_to_tsv: pathlib.Path) -> pd.DataFrame:
    return pd.read_table(path_to_tsv)


def import_chrom_sizes(path_to_chrom_sizes: Union[pathlib.Path, None]) -> Union[pd.DataFrame, None]:
    if path_to_chrom_sizes is None:
        return None

    return pd.read_table(path_to_chrom_sizes, names=["chrom", "length"])


def compute_weighted_avg_and_std(data: pd.Series, weights: Union[pd.Series, None]) -> Tuple[float, float]:
    if weights is not None:
        assert len(data) == len(weights)

    avg = np.average(data, weights=weights)
    variance = np.average((data - avg) ** 2, weights=weights)

    return avg, np.sqrt(variance)


def aggregate_scc(df: pd.DataFrame, weights: Union[pd.DataFrame, None]) -> Tuple[np.ndarray, np.ndarray]:
    conditions = df["cond1"].unique().tolist()

    scc = np.zeros([len(conditions), len(conditions)])
    std = np.zeros([len(conditions), len(conditions)])
    for i in range(len(conditions)):
        for j in range(i, len(conditions)):
            if i == j:
                scc[i, j] = 1.0
                continue

            df1 = df[(df["cond1"] == conditions[i]) & (df["cond2"] == conditions[j])]
            corr, stddev = compute_weighted_avg_and_std(df1["scc"], weights)

            scc[i, j] = corr
            scc[j, i] = corr

            std[i, j] = stddev
            std[j, i] = stddev

    return scc, std


def plot_heatmap(scc: np.ndarray, std: np.ndarray, conditions: List[str], title: str) -> plt.Figure:
    fig, ax = plt.subplots(1, 1)

    min_val = scc.min()
    delta = scc.max() - min_val

    ax.imshow(scc, cmap="Reds")
    for (j, i), label in np.ndenumerate(scc):
        if label - min_val >= (delta * 0.85):
            color = "white"
        else:
            color = "black"

        if i == j:
            label = f"{label:.3f}"
        else:
            label = f"{label:.3f}\n±\n{std[i, j]:.3f}"
        ax.text(i, j, label, ha="center", va="center", color=color)

    ax.set_xticks(range(len(conditions)))
    ax.set_yticks(range(len(conditions)))

    ax.set_xticklabels(conditions, rotation=30)
    ax.set_yticklabels(conditions)
    ax.set(title=title)

    plt.tight_layout()

    return fig


def main():
    args = vars(make_cli().parse_args())

    output_prefix = args["output_prefix"]
    output_prefix.parent.mkdir(exist_ok=True)

    df = import_tsv(args["tsv"])
    chroms = df["chrom"].unique().tolist()
    chrom_sizes = import_chrom_sizes(args.get("chrom-sizes"))

    if chrom_sizes is None:
        weights = None
    else:
        weights = chrom_sizes.T[[chroms]].T

    scc, std = aggregate_scc(df, weights)

    conditions = df["cond1"].unique().tolist()
    fig = plot_heatmap(scc, std, conditions, "HiC replicate SCC (weighted average)")

    fig.savefig(output_prefix.with_suffix(".svg"))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600)


if __name__ == "__main__":
    main()
