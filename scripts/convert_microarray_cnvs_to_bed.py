#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import bioframe as bf
import numpy as np
import pandas as pd
import subprocess as sp
import tempfile
import shutil
import logging


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--probe-ids",
        nargs="+",
        required=True,
        type=pathlib.Path,
        help="File mapping SNP probe IDs to genomic coordinates (e.g. GPL3718 and GPL3720).",
    )

    cli.add_argument(
        "--liftover-chain", type=pathlib.Path, help="Path to the chain file to use when lifting-over coordinates."
    )

    cli.add_argument(
        "--path-to-liftover", type=pathlib.Path, default=pathlib.Path("liftOver"), help="Path to UCSC liftOver."
    )

    cli.add_argument("--fill-gaps", action="store_true", default=False, help="Fill gaps between SNP probes.")

    cli.add_argument(
        "--chrom-sizes",
        type=pathlib.Path,
        required="--fill-gaps" in sys.argv,
        help="Path to the .chrom.sizes file to use when filling gaps between SNPs",
    )

    cli.add_argument("tsv", type=pathlib.Path, help="TSV with CNVs (e.g. *_copynumber.txt.gz from GSE19920).")

    return cli


def import_probes(path_to_txt: pathlib.Path) -> pd.DataFrame:
    df = pd.read_table(path_to_txt, comment="#", low_memory=False)
    df = df[~df["ID"].str.startswith("AFFX-")]  # Drop control probes
    df["Chromosome"] = "chr" + df["Chromosome"].str.removeprefix("chr")
    return df


def extract_probe_lengths(df: pd.DataFrame) -> np.ndarray:
    # Entries in the SEQUENCE column look something like this: ggatattgtgtgagga[A/G]taagcccacctgtggt
    probes = df["SEQUENCE"].str.replace(r"[\[\]A-Z\/]+", "", n=1, regex=True)
    lengths = probes.str.len() + 1
    return np.array(lengths, dtype=int)


def probes_to_bed(df: pd.DataFrame) -> pd.DataFrame:
    probe_lengths = extract_probe_lengths(df)

    df = df[["Chromosome", "Physical Position", "ID", "Strand"]].copy()
    df.columns = ["chrom", "start", "name", "strand"]

    df["start"] = df["start"].astype(int)
    df["end"] = df["start"] + probe_lengths
    df["score"] = "."

    return bf.sort_bedframe(df[bf.SCHEMAS["bed6"]])


def run_liftover(bed: pd.DataFrame, path_to_chain: pathlib.Path, path_to_liftover: pathlib.Path) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile(mode="rt") as hits, tempfile.NamedTemporaryFile(mode="rt") as miss:
        cmd = [path_to_liftover, "stdin", str(path_to_chain), str(hits.name), str(miss.name)]

        logging.info("Running liftOver...")
        with sp.Popen(cmd, stdin=sp.PIPE, stderr=sp.PIPE) as liftover:
            bed.to_csv(liftover.stdin, sep="\t", header=False, index=False)
            liftover.communicate(timeout=100)

            if (code := liftover.returncode) != 0:
                print(liftover.stderr, file=sys.stderr)
                raise RuntimeError(f"liftOver terminated with code {code}")

        logging.info("liftOver terminated without errors!")

        for misses, _ in enumerate(miss.file):
            pass

        if misses != 0:
            logging.warning("liftOver failed to map %s entries", misses)

        return pd.read_table(hits.file, names=bed.columns)


def find_liftover_or_fail(path: pathlib.Path) -> str:
    liftover = shutil.which(path)
    if liftover is None:
        raise RuntimeError(f"which: unable to find {path} in your PATH")

    return liftover


def fill_gaps(df: pd.DataFrame, chrom_sizes: dict) -> pd.DataFrame:
    # Loop over chromosomes, and fill gaps between consecutive intervals as shown below:
    # Given two adjacent intervals in the sorted dataframe df:
    #  chr1 10 20 ...
    #  chr1 50 60 ...
    # Fill the gap like so:
    #  chr1 10   30.5
    #  chr1 30.5 60
    dfs = []
    for chrom, df2 in df.groupby("chrom"):
        if chrom not in chrom_sizes:
            continue
        df1 = df2.shift(1)  # prev
        df3 = df2.shift(-1)  # next

        new_starts = (df1["end"] + df2["start"]).to_numpy() / 2
        new_ends = (df2["end"] + df3["start"]).to_numpy() / 2

        # Cap chromosome
        new_starts[0] = 0
        new_ends[-1] = chrom_sizes[chrom]

        df2["start"] = np.round(new_starts).astype(int)
        df2["end"] = np.round(new_ends).astype(int)
        dfs.append(df2)

    cols = [col for col in df.columns if col not in bf.SCHEMAS["bed6"]]
    df1 = bf.sort_bedframe(bf.merge(pd.concat(dfs), min_dist=0, on=cols)).drop(columns=["n_intervals"])

    # Recover probe IDs and store them as a ";" separated string in the name field
    df2 = bf.overlap(df1, df[bf.SCHEMAS["bed4"]])[["chrom", "start", "end", "name_"]]
    probe_ids = df2.groupby(by=bf.SCHEMAS["bed3"])["name_"].transform(lambda pid: ";".join(pid)).drop_duplicates()

    df1["name"] = probe_ids.reset_index(drop=True).tolist()

    cols = bf.SCHEMAS["bed4"] + [col for col in df1.columns if col not in bf.SCHEMAS["bed4"]]
    return df1[cols]


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    if args.get("liftover_chain"):
        args["path_to_liftover"] = find_liftover_or_fail(args["path_to_liftover"])

    probe_db = pd.concat([import_probes(probe_file) for probe_file in args["probe_ids"]])
    logging.info(
        "Imported %s unique probes from file(s) %s", len(probe_db), ", ".join((str(p) for p in args["probe_ids"]))
    )
    probes = probes_to_bed(probe_db)

    if chain := args.get("liftover_chain"):
        probes = run_liftover(probes, chain, args["path_to_liftover"])

    cnvs = pd.read_table(args["tsv"])
    cnvs = bf.sort_bedframe(
        cnvs.merge(probes, left_on="ProbeSet", right_on="name").drop(columns=["ProbeSet", "Chromosome", "Position"])
    )

    if num_nans := cnvs.isnull().values.sum() != 0:
        raise RuntimeError(f"Failed to convert {num_nans} records")

    if args["fill_gaps"]:
        logging.info("Extending itervals to fill gaps betwen SNP probes...")
        chrom_sizes = pd.read_table(args["chrom_sizes"], names=["chrom", "size"]).set_index("chrom")["size"].to_dict()
        cnvs = fill_gaps(cnvs, chrom_sizes)

    cnvs.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    setup_logger()
    main()
