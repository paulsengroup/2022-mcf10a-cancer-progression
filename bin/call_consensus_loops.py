#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import subprocess as sp
import sys
import tempfile
import warnings

import bioframe as bf
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedpe",
        nargs="+",
        type=pathlib.Path,
        help="Path to two or more BEDPE file(s) with the list of loops to be processed.",
    )

    cli.add_argument("--loop-resolution", required=True, type=int, help="Lowest resolution used to call loops.")

    cli.add_argument(
        "--min-overlap",
        type=float,
        default=0.01,
        help="Minimum relative overlap required in order for two dots to be considered the same.",
    )

    cli.add_argument("--vote-type", default="majority", type=str, choices={"majority", "all"})

    cli.add_argument(
        "--skip-header",
        action="store_true",
        default=False,
        help="Skip file header.",
    )

    return cli


def run_pair_to_pair(df1, df2, type="both", padding=0, frac=0.01):
    with tempfile.NamedTemporaryFile("wt") as bedpe1, tempfile.NamedTemporaryFile("wt") as bedpe2:
        df1.to_csv(bedpe1.name, sep="\t", header=False, index=False)
        df2.to_csv(bedpe2.name, sep="\t", header=False, index=False)

        cmd = [
            "pairToPair",
            "-a",
            bedpe1.name,
            "-b",
            bedpe2.name,
            "-is",
            "-type",
            type,
            "-f",
            str(frac),
        ]
        if padding != 0:
            cmd.extend(("-slop", str(padding)))
        with sp.Popen(cmd, stdout=sp.PIPE) as bedtools:
            df = pd.read_table(bedtools.stdout, usecols=list(range(6)), names=bf.SCHEMAS["bedpe"][:6])
            bedtools.communicate(timeout=20)

        if (code := bedtools.returncode) != 0:
            print(bedtools.stderr, file=sys.stderr)
            raise RuntimeError(f"bedtools terminated with code {code}")

        return df.drop_duplicates(keep="first")


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


def estimate_loop_resolution(df: pd.DataFrame, num_samples: int = 100) -> int:
    df1 = df.sample(num_samples, random_state=len(df))
    return (df1["end1"] - df1["start1"]).median()


def main():
    args = vars(make_cli().parse_args())

    assert len(args["bedpe"]) > 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dfs = {}
        for p in args["bedpe"]:
            df = bf.read_table(p, schema=bf.SCHEMAS["bedpe"][:6], skiprows=int(args["skip_header"]))
            loop_size = estimate_loop_resolution(df)
            assert loop_size not in dfs
            dfs[loop_size] = pad_loops(df, args["loop_resolution"])

    # Sort loops by resolution
    dfs = [v for k, v in reversed(sorted(dfs.items()))]

    # Pair loops
    df1 = dfs[0]
    dfs_ = []
    for df2 in dfs[1:]:
        dfs_.append(run_pair_to_pair(df1, df2, frac=args["min_overlap"]))

    # Count votes
    df = pd.concat(dfs_).groupby(bf.SCHEMAS["bedpe"][:6]).size() + 1
    df = df.rename("votes").to_frame().reset_index()

    if args["vote_type"] == "majority":
        df = df[df["votes"] >= len(dfs_) / 2]
    else:
        assert args["vote_type"] == "all"
        df = df[df["votes"] == len(dfs_)]

    df.to_csv(sys.stdout, header=False, index=False, sep="\t")


if __name__ == "__main__":
    main()
