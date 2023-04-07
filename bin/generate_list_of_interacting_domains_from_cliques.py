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
        description="Convert the output of robomics/call_tad_cliques to a BEDPE with pairs of interacting domains."
    )

    cli.add_argument("domains", type=existing_file, help="Path to a BED file with the list of domains.")
    cli.add_argument("cliques", type=existing_file, help="Path to a TSV file with the list of cliques.")

    return cli


def read_domains(domains: pathlib.Path) -> pd.DataFrame:
    return pd.read_table(domains, names=["chrom", "start", "end", "id"]).set_index("id")


def read_cliques(cliques: pathlib.Path) -> pd.DataFrame:
    df = pd.read_table(cliques).rename(columns={"name": "clique"})
    assert df.columns.tolist() == ["clique", "tad_ids", "size"]

    df["tad_ids"] = df["tad_ids"].str.split(",")
    return df.set_index("clique")


def extract_domain_pairs_from_cliques(cliques: pd.DataFrame) -> pd.DataFrame:
    """
    Generate a dataframe with the list of pairs of interacting domains.
    The list does not contain duplicates.
    """
    domain_pairs = set()

    cliques["tad_ids"].apply(
        lambda doms: domain_pairs.update(
            tuple(map(int, pair)) for pair in itertools.combinations_with_replacement(doms, 2)
        )
    )
    assert len(domain_pairs) != 0

    return pd.DataFrame(list(domain_pairs), columns=["tad1_id", "tad2_id"])


def map_domain_coords_to_domain_pairs(domains: pd.DataFrame, coords: pd.DataFrame) -> pd.DataFrame:
    df = domains.merge(coords, left_on="tad1_id", right_index=True, suffixes=("", "1"))
    df = df.merge(coords, left_on="tad2_id", right_index=True, suffixes=("1", "2"))

    return df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].sort_values(
        by=["chrom1", "start1", "chrom2", "start2"]
    )


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

    domain_pairs = extract_domain_pairs_from_cliques(cliques)
    domain_pairs = map_domain_coords_to_domain_pairs(domain_pairs, domains)

    domain_pairs.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
