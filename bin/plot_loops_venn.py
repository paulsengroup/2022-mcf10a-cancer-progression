#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import warnings

import bioframe as bf
import matplotlib.pyplot as plt
import matplotlib_venn
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli.add_argument(
        "loops",
        nargs="+",
        type=pathlib.Path,
        help="Path to three BEDPE with the list of loops to plot.",
    )

    cli.add_argument("lowest-resolution", type=positive_int, help="Lowest resolution used for loop calling.")

    cli.add_argument(
        "--labels",
        type=str,
        nargs="+",
        help="List of labels to use instead of file names.",
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
        help="Overwrite existing files (if any).",
    )

    return cli


def import_bedpe(path_to_bedpe: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cols = bf.SCHEMAS["bedpe"][:6]
        return bf.read_table(path_to_bedpe, schema="bedpe")[cols]


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


def pad_loops(df: pd.DataFrame, target_size: int) -> pd.DataFrame:
    df = df.copy()
    df["span"] = df["end1"] - df["start1"]
    dfs = []

    for span, grp in df.groupby("span"):
        if target_size < span:
            dfs.append(grp)
            continue

        padding = target_size - span
        grp = bf.expand(grp, padding // 2, cols=("chrom1", "start1", "end1"))
        grp = bf.expand(grp, padding // 2, cols=("chrom2", "start2", "end2"))

        dfs.append(grp.drop(columns="span"))

    return pd.concat(dfs)


def overlap_bedpe(*dfs: pd.DataFrame) -> pd.DataFrame:
    df1 = dfs[0]

    for df2 in dfs[1:]:
        # Overlap first triplet
        df3 = bf.overlap(df1, df2, cols1=("chrom1", "start1", "end1"), cols2=("chrom1", "start1", "end1")).dropna()

        data = []

        # Overlap second triplet
        for key, grp in df3.groupby(["chrom1", "start1", "end1"]):
            chrom1, start1, end1 = key
            df4 = grp[["chrom2", "start2", "end2"]].copy()
            df5 = grp[["chrom2_", "start2_", "end2_"]].copy()

            df4["start2"] = df4["start2"].astype(int)
            df5["start2_"] = df5["start2_"].astype(int)

            df4["end2"] = df4["end2"].astype(int)
            df5["end2_"] = df5["end2_"].astype(int)

            overlap = bf.closest(df4, df5, cols1=("chrom2", "start2", "end2"), cols2=("chrom2_", "start2_", "end2_"))

            if len(overlap) == 0:
                continue

            data.append(
                [chrom1, start1, end1, overlap.loc[0, "chrom2"], overlap.loc[0, "start2"], overlap.loc[0, "end2"]]
            )

        df1 = pd.DataFrame.from_records(data, columns=bf.SCHEMAS["bedpe"][:6])
    return df1


def main():
    args = vars(make_cli().parse_args())

    labels = args.get("labels")
    paths_to_loops = args["loops"]
    if len(paths_to_loops) != 3:
        raise RuntimeError(f"Expected path to 3 loop annotations, found {len(paths_to_loops)}")

    if labels is None:
        labels = [p.name for p in paths_to_loops]
    elif len(labels) != len(paths_to_loops):
        raise RuntimeError(
            f"Mismatch in the number of files and labels: expected {len(paths_to_loops)} labels, found {len(labels)}"
        )

    loops = {k: pad_loops(import_bedpe(path), args["lowest-resolution"]) for k, path in zip(labels, paths_to_loops)}

    a = loops[labels[0]]
    b = loops[labels[1]]
    c = loops[labels[2]]

    a_vs_b = overlap_bedpe(a, b)
    a_vs_c = overlap_bedpe(a, c)
    b_vs_c = overlap_bedpe(b, c)
    shared = overlap_bedpe(a, b, c)

    ab = len(a_vs_b) - len(shared)
    bc = len(b_vs_c) - len(shared)
    ac = len(a_vs_c) - len(shared)

    a_only = len(a) - ab - ac - len(shared)
    b_only = len(b) - ab - bc - len(shared)
    c_only = len(c) - ac - bc - len(shared)

    subsets = {"100": a_only, "010": b_only, "001": c_only, "110": ab, "101": ac, "011": bc, "111": len(shared)}

    fig, ax = plt.subplots(1, 1)
    matplotlib_venn.venn3(subsets, set_labels=labels, ax=ax)

    save_plot_to_file(fig, args["output_prefix"], args["force"])


if __name__ == "__main__":
    main()
