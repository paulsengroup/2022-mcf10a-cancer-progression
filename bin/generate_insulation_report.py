#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings
from typing import Dict, Tuple

import bioframe as bf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "insulation",
        type=pathlib.Path,
        help="Path to the aggregated insulation score bedgraph.",
    )

    cli.add_argument(
        "domains",
        type=pathlib.Path,
        help="Path to a BED file with the list of TADs to be used in the report (usually WT).",
    )

    cli.add_argument(
        "--blacklist", type=pathlib.Path, help="Path to a BED file with the list of regions to be excluded."
    )

    cli.add_argument(
        "--blacklist-padding", type=int, default=100_000, help="Number of bps to add as padding to blacklisted regions."
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        help="Output prefix.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def save_plot_to_file(fig: plt.Figure, outprefix: pathlib.Path, force: bool, close_after_save: bool = True) -> None:
    png = outprefix.with_suffix(".png")
    svg = outprefix.with_suffix(".svg")
    if not force:
        handle_path_collisions(png, svg)

    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight")
    if close_after_save:
        plt.close(fig)


def import_bedgraph(path_to_bedgraph: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path_to_bedgraph, schema="bedGraph").set_index(["chrom", "start", "end"])


def import_insulation_scores(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_table(path).set_index(["chrom", "start", "end"])


def import_tads(path: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path, schema="bed3").set_index(["chrom", "start", "end"])


def detect_path_collisions(*args: pathlib.Path) -> None:
    for path in args:
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite file {path}. Pass --force to overwrite existing file(s).")


def compute_insulation_std(scores: pd.DataFrame, mappings: Dict) -> pd.DataFrame:
    stds = []
    for _, row in scores.iterrows():
        std = 0.0
        for columns in mappings.values():
            std += row.loc[columns].std(skipna=False)
        stds.append(std)

    scores["std"] = stds
    return scores.sort_values("std", ascending=False)


def map_insulation_to_tads(domains: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    df1 = domains.reset_index().copy()

    df1["end"] = df1["start"] + 1

    df1 = (
        bf.overlap(df1, scores.reset_index(), suffixes=("", "_"))
        .drop(columns=["chrom_", "start_", "end_"])
        .rename(columns=lambda c: c.rstrip("_"))
    )

    df2 = (
        bf.overlap(domains.reset_index(), df1, suffixes=("", "_"))
        .drop(columns=["chrom_", "start_", "end_"])
        .rename(columns=lambda c: c.rstrip("_"))
        .set_index(["chrom", "start", "end"])
    )

    return df2


def compute_axis_scale(scores: pd.DataFrame) -> Tuple[float, float]:
    lb = np.inf
    ub = -np.inf
    num_cols = len(scores.columns)
    for i in range(num_cols):
        for j in range(i + 1, num_cols):
            cond1, cond2 = scores.columns[i], scores.columns[j]
            df = scores[[cond1, cond2]].dropna()
            lb = min(lb, df[cond1].min(), df[cond2].min())
            ub = max(ub, df[cond1].max(), df[cond2].max())

    ub += abs(ub) * 0.05
    lb -= abs(lb) * 0.05

    return lb, ub


def plot_scatters(fig, axs, lb, ub, scores: pd.DataFrame):
    num_cols = len(scores.columns)
    for i in range(num_cols):
        for j in range(i + 1, num_cols):
            cond1, cond2 = scores.columns[i], scores.columns[j]
            df = scores.copy()
            df["delta"] = (df[cond1] - df[cond2]).abs()
            axs[i][j].scatter(df[cond1], df[cond2], color="blue", alpha=0.1)
    for i in range(num_cols):
        for j in range(num_cols):
            ax = axs[i][j]

            if i == 0:
                ax.set(title=scores.columns[j])

            if j == i + 1:
                ax.set(ylabel=scores.columns[i])

            if i >= j:
                fig.delaxes(ax)
            else:
                ax.set(xlim=(lb, ub), ylim=(lb, ub))


def main():
    args = vars(make_cli().parse_args())

    paths_to_scores = args["insulation"]
    paths_to_domains = args["domains"]

    scores = import_insulation_scores(paths_to_scores)
    domains = import_tads(paths_to_domains)

    if args["blacklist"] is not None:
        df1 = domains.copy()

        for col in df1.index.names:
            df1[col] = df1.index.get_level_values(col)

        df1 = bf.expand(df1, args["blacklist_padding"])
        domains = bf.setdiff(df1, bf.read_table(args["blacklist"], schema="bed")).drop(
            columns=["chrom", "start", "end"]
        )

    fig, ax = plt.subplots(1, 1)
    ax.hist(scores["std"], bins=25, log=True)
    ax.set(title="Insulation score STD", xlabel="std", ylabel="frequency")
    outprefix = pathlib.Path(str(args["output_prefix"]) + "_hist")
    plt.tight_layout()
    save_plot_to_file(fig, outprefix, args["force"])

    scores = map_insulation_to_tads(domains, scores)

    grid_size = len(scores.columns) - 1
    fig, axs = plt.subplots(grid_size, grid_size, figsize=(grid_size * 6.4, grid_size * 6.4))
    lb, ub = compute_axis_scale(scores)
    plot_scatters(fig, axs, lb, ub, scores.drop(columns=["std"]))

    plt.tight_layout()

    outprefix = pathlib.Path(str(args["output_prefix"]) + "_scatter")
    save_plot_to_file(fig, outprefix, args["force"])


if __name__ == "__main__":
    main()
