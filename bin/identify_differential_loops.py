#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import itertools
import pathlib
import sys
import warnings

import bioframe as bf
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
        help="Path to two or more BEDPE with the list of loops to compare.",
    )

    cli.add_argument(
        "lowest-resolution",
        type=positive_int,
        help="Lowest resolution used for loop calling.",
    )

    cli.add_argument(
        "--labels",
        type=str,
        nargs="+",
        help="List of labels to use instead of file names.",
    )

    return cli


def import_bedpe(path_to_bedpe: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cols = bf.SCHEMAS["bedpe"][:6]
        return bf.read_table(path_to_bedpe, schema="bedpe")[cols]


def detect_path_collisions(*args: pathlib.Path) -> None:
    for path in args:
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite file {path}. Pass --force to overwrite existing file(s).")


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
        df3 = bf.overlap(
            df1,
            df2,
            cols1=("chrom1", "start1", "end1"),
            cols2=("chrom1", "start1", "end1"),
        ).dropna()

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

            overlap = bf.closest(
                df4,
                df5,
                cols1=("chrom2", "start2", "end2"),
                cols2=("chrom2_", "start2_", "end2_"),
            )

            if len(overlap) == 0:
                continue

            data.append(
                [
                    chrom1,
                    start1,
                    end1,
                    overlap.loc[0, "chrom2"],
                    overlap.loc[0, "start2"],
                    overlap.loc[0, "end2"],
                ]
            )

        df1 = pd.DataFrame.from_records(data, columns=bf.SCHEMAS["bedpe"][:6])
    return df1


def main():
    args = vars(make_cli().parse_args())

    labels = args.get("labels")
    paths_to_loops = args["loops"]
    if labels is None:
        labels = [p.name for p in paths_to_loops]
    elif len(labels) != len(paths_to_loops):
        raise RuntimeError(
            f"Mismatch in the number of files and labels: expected {len(paths_to_loops)} labels, found {len(labels)}"
        )

    loops = {k: pad_loops(import_bedpe(path), args["lowest-resolution"]) for k, path in zip(labels, paths_to_loops)}

    print_header = True
    for (k1, l1), (k2, l2) in itertools.product(loops.items(), repeat=2):
        if k1 == k2:
            continue

        df = overlap_bedpe(l1, l2)
        df = df[~df["has_overlap"]].drop(columns="has_overlap")
        df["cond1"] = k1
        df["cond2"] = k2
        df.to_csv(sys.stdout, sep="\t", header=print_header, index=True)
        print_header = False


if __name__ == "__main__":
    main()
