#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys
import warnings
from typing import Dict, Tuple, Union

import bioframe as bf
import cooler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "coolers",
        type=pathlib.Path,
        nargs="+",
        help="Path to two or more cooler files (URI syntax supported).",
    )

    cli.add_argument(
        "domains",
        type=pathlib.Path,
        help="Path to a BED file with the list of TADs to be used for plotting (usually WT).",
    )

    cli.add_argument("--diagonals-to-mask", type=int, default=5, help="Number of diagonal to mask.")

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        help="Output prefix.",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma separated list of labels to use instead of file names.",
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


def import_tads(path: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path, schema="bed3").set_index(["chrom", "start", "end"])


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


def get_interactions(coords, sel: cooler.api.RangeSelector2D, diagonals_to_mask) -> Union[int, float]:
    chrom, start, end = coords
    m = sel.fetch(f"{chrom}:{start}-{end}").tolil()
    for i in range(min(m.shape[0], diagonals_to_mask)):
        m.setdiag(0, i)

    m = m.tocoo()
    return m.sum() / m.max()


def main():
    args = vars(make_cli().parse_args())

    paths_to_domains = args["domains"]
    cooler_uris = args["coolers"]
    labels = args.get("labels")

    if labels is None:
        labels = [cooler.Cooler(str(p)).filename for p in cooler_uris]
    else:
        labels = labels.split(",")
        if len(labels) != len(cooler_uris):
            raise RuntimeError(
                f"Mismatch in the number of files and labels: expected {len(cooler_uris)} labels, found {len(labels)}"
            )

    domains = import_tads(paths_to_domains)

    for key, uri in zip(labels, args["coolers"]):
        sel = cooler.Cooler(str(uri)).matrix(sparse=True, balance=False)
        domains[key] = domains.index.to_frame().apply(
            get_interactions,
            args=(sel, args["diagonals_to_mask"]),
            axis="columns",
        )
    grid_size = len(domains.columns)
    fig, axs = plt.subplots(grid_size, grid_size, figsize=(grid_size * 6.4, grid_size * 6.4))
    lb, ub = compute_axis_scale(domains)
    plot_scatters(fig, axs, lb, ub, domains)

    plt.tight_layout()

    outprefix = pathlib.Path(str(args["output_prefix"]) + "_scatter")
    save_plot_to_file(fig, outprefix, args["force"])


if __name__ == "__main__":
    main()