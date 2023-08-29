#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import warnings
from typing import Dict, List

import bioframe as bf
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "domains",
        nargs="+",
        type=pathlib.Path,
        help="Path to two or more BED files to overlap.",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        type=pathlib.Path,
        help="Output prefix.",
    )

    cli.add_argument("--labels", type=str, help="Comma separated list of labels to use instead of file names.")

    cli.add_argument(
        "--blacklist", type=pathlib.Path, help="Path to a BED file with the list of regions to be excluded."
    )
    cli.add_argument("--title", type=str, help="Plot title.")

    cli.add_argument(
        "--blacklist-padding", type=int, default=0, help="Number of base-pairs used to pad blacklisted intervals."
    )
    cli.add_argument("--overlap-threshold", type=float, default=0.8)

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def import_bed(path_to_bed: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path_to_bed, schema="bed3").set_index(["chrom", "start", "end"])


def import_tads(paths: List[pathlib.Path], labels: List[str]) -> Dict[str, pd.DataFrame]:
    assert len(paths) == len(labels)
    dfs = {}
    for label, path in zip(labels, paths):
        dfs[label] = import_bed(path)

    return dfs


def compute_relative_overlap(df1, df2):
    df = bf.overlap(df1, df2, return_overlap=True)
    df["overlap"] = df["overlap_end"] - df["overlap_start"]
    df["rel_overlap"] = df["overlap"] / (np.maximum(df["end"], df["end_"]) - np.minimum(df["start"], df["start_"]))
    df = df.sort_values(["chrom", "start", "end", "rel_overlap"]).groupby(["chrom", "start", "end"]).last()

    df["overlap"] = df["overlap"].fillna(0)
    df["rel_overlap"] = df["rel_overlap"].fillna(0)
    return df


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


def main():
    args = vars(make_cli().parse_args())

    labels = args.get("labels")
    paths_to_domains = args["domains"]

    if labels is None:
        labels = [p.name for p in paths_to_domains]
    else:
        labels = labels.split(",")
        if len(labels) != len(paths_to_domains):
            raise RuntimeError(
                f"Mismatch in the number of files and labels: expected {len(paths_to_domains)} labels, found {len(labels)}"
            )

    domains = import_tads(paths_to_domains, labels)

    if args["blacklist"] is not None:
        blacklist = bf.read_table(args["blacklist"], schema="bed")
        for key, df in domains.items():
            df = df.copy()

            for col in df.index.names:
                df[col] = df.index.get_level_values(col)

            df = bf.expand(df, args["blacklist_padding"])
            domains[key] = bf.setdiff(df, blacklist).drop(columns=["chrom", "start", "end"])

    fig, ax = plt.subplots(1, 1)

    labels = []
    high_overlap_threshold = args["overlap_threshold"]
    for i, (cond1, df1) in enumerate(domains.items()):
        for cond2, df2 in list(domains.items())[i + 1 :]:
            df = compute_relative_overlap(df1.reset_index(), df2.reset_index())
            ax.hist(df["rel_overlap"], bins=50, label=f"{cond1}_{cond2}", histtype="step")

            high_overlap = (df["rel_overlap"] >= high_overlap_threshold).sum() / len(df)
            labels.append(f"{cond1}_{cond2}: {high_overlap:.2f}")
    color_palette = sns.color_palette().as_hex()
    artists = [Rectangle((0, 0), 3, 1, facecolor=color_palette[i], edgecolor="black") for i, _ in enumerate(labels)]

    ax.legend(
        handles=artists,
        labels=labels,
        loc="upper left",
        title=f"TADs >= {high_overlap_threshold}% overlap",
    )

    if args["title"] is None:
        args["title"] = f"TAD overlap"
    fig.suptitle(args["title"])

    plt.tight_layout()

    save_plot_to_file(fig, args["output_prefix"], args["force"])


if __name__ == "__main__":
    mpl.rcParams.update(
        {
            "axes.titlesize": 10,
            "axes.labelsize": 22,
            "legend.fontsize": 17,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        }
    )
    main()
