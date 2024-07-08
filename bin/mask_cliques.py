#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import bioframe as bf
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=existing_file,
        help="Path to the TSV file with the list of cliques.",
    )
    cli.add_argument(
        "bed",
        type=existing_file,
        help="Path to the BED file with the list of domains.",
    )

    cli.add_argument(
        "--mask",
        type=existing_file,
        help="Path to a BED3+ file with the list of regions to be masked.",
    )

    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    cliques = pd.read_table(args["tsv"])
    if args["mask"] is None:
        cliques.to_csv(sys.stdout, sep="\t", header=True, index=False)
        sys.exit(0)

    domains = pd.read_table(args["bed"], names=["chrom", "start", "end", "id"])
    mask = pd.read_table(args["mask"], names=["chrom", "start", "end"], usecols=list(range(3)))

    domains = bf.setdiff(domains, mask)
    domains["id"] = domains["id"].astype(str)

    cliques["tad_ids"] = cliques["tad_ids"].str.split(",")
    cliques = cliques.explode("tad_ids")

    cliques = cliques.merge(domains, left_on="tad_ids", right_on="id").drop(columns="id")

    rows = []
    for name, df in cliques.groupby("name"):
        if df.isnull().any().any():
            continue
        rows.append([name, ",".join(df["tad_ids"]), len(df)])

    df = pd.DataFrame(rows, columns=["name", "tad_ids", "size"])
    df.to_csv(sys.stdout, sep="\t", header=True, index=False)
