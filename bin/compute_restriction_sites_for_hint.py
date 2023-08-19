#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse

import bioframe as bf
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("fna", type=str)
    cli.add_argument("enzymes", nargs="+", type=str)

    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    dfs = []

    fna = bf.load_fasta(args["fna"])
    for enz in args["enzymes"]:
        dfs.append(bf.digest(fna, enz))

    df = pd.concat(dfs)[["chrom", "start"]].drop_duplicates(keep="first")

    groups = df.groupby("chrom")["start"].aggregate(list)
    for chrom, sites in groups.apply(lambda l: " ".join(str(x) for x in sorted(l))).items():
        print(chrom, sites)
