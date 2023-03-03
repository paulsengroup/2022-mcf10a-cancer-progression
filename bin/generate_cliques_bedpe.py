#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT


import argparse
import itertools
import pathlib
import sys

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser(
        description="Convert the output of robomics/call_tad_cliques to a BEDPE with pairs of domains taking part in maximal cliques."
    )

    cli.add_argument("domains", type=existing_file, help="Path to a BED file with the list of domains.")
    cli.add_argument("cliques", type=existing_file, help="Path to a TSV file with the list of cliques.")

    return cli


def read_domains(domains) -> pd.DataFrame:
    return pd.read_table(domains, names=["chrom", "start", "end", "id"]).set_index("id")


def read_cliques(cliques) -> pd.DataFrame:
    df = pd.read_table(cliques).rename(columns={"name": "clique"})
    assert df.columns.tolist() == ["clique", "tad_ids", "size"]
    df = df.set_index("clique")

    domains = df["tad_ids"].str.split(",")

    pairs = domains.apply(lambda doms: tuple(itertools.product(doms, repeat=2)))

    data = []
    for idx, row in pairs.items():
        for pair in row:
            data.append(tuple([idx, *pair]))

    df = pd.DataFrame(data, columns=["clique", "tad1_id", "tad2_id"]).infer_objects()
    df["tad1_id"] = pd.to_numeric(df["tad1_id"])
    df["tad2_id"] = pd.to_numeric(df["tad2_id"])

    return df


def main():
    args = vars(make_cli().parse_args())

    domains = read_domains(args["domains"])
    cliques = read_cliques(args["cliques"])

    if len(domains) == 0:
        src = args["domains"]
        raise RuntimeError(f"Unable to read any record from {src}")

    if len(cliques) == 0:
        src = args["cliques"]
        raise RuntimeError(f"Unable to read any record from {src}")

    df = cliques.merge(domains, left_on="tad1_id", right_index=True, suffixes=("", "1"))
    df = df.merge(domains, left_on="tad2_id", right_index=True, suffixes=("1", "2"))

    df = df[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "clique"]].sort_values(
        by=["chrom1", "start1", "chrom2", "start2"]
    )
    df = df.sort_values(by=["clique"], kind="stable", key=lambda idx: idx.str.removeprefix("CLIQUE_#").astype(int))

    df.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()