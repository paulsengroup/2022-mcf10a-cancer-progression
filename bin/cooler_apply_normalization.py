#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import multiprocessing as mp
import pathlib
import shutil
import subprocess as sp
import sys
from typing import Union

import cooler
import pandas as pd
from cooler.core import CSRReader, DirectRangeQuery2D


def parse_cooler_uri(uri):
    cf = cooler.Cooler(uri)
    return cf.filename, cf.root.lstrip("/")


def printable_chrom(chrom):
    if chrom is None:
        return "all"
    return str(chrom)


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "cooler",
        type=pathlib.Path,
        help="Path to a .cool file (URI syntax supported).",
    )
    cli.add_argument(
        "--norm-name",
        type=str,
        default="weight",
        help="Name of the normalization to be applied.",
    )
    cli.add_argument(
        "--weight-type",
        type=str,
        default="auto",
        choices={"auto", "divisive", "multiplicative"},
        help="Specify whether weights are multiplicative or divisive.\n"
        "By default, weight type is inferred from the dataset name.",
    )
    cli.add_argument(
        "--output-name",
        type=pathlib.Path,
        required=True,
        help="Output file name.\n"
        "When output name has extension other than .cool,\n"
        "normalized interactions are stored in BEDPE format.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )
    return cli


def weights_are_divisive(weight_name: str) -> bool:
    divisive_weight_names = {
        "VC",
        "INTER_VC",
        "GW_VC",
        "VC_SQRT",
        "KR",
        "INTER_KR",
        "GW_KR",
        "SCALE",
        "INTER_SCALE",
        "GW_SCALE",
    }

    return weight_name in divisive_weight_names


def make_annotator(bins, weight_dset, join, divisive):
    # Adapted from: https://github.com/open2c/cooler/blob/master/cooler/cli/dump.py
    def annotator(chunk):
        df = cooler.annotate(chunk, bins[[weight_dset]])
        w1 = df[f"{weight_dset}1"]
        w2 = df[f"{weight_dset}2"]

        if divisive:
            w1 = 1 / w1
            w2 = 1 / w2

        chunk["balanced"] = w1 * w2 * chunk["count"]

        if join:
            chunk = cooler.annotate(chunk, bins[["chrom", "start", "end"]], replace=True)

        return chunk

    return annotator


def get_compressor(compression_level: int, nproc: int) -> Union[list, None]:
    if cmd := shutil.which("pigz"):
        return [cmd, f"-{compression_level}", "-p", str(min(16, nproc))]
    if cmd := shutil.which("gzip"):
        return [cmd, f"-{compression_level}"]
    return None


def to_cooler(
    chunks,
    bins: pd.DataFrame,
    output_cooler_uri: pathlib.Path,
    assembly: str = "unknown",
):
    chunks = (chunk.rename(columns={"balanced": "count"}) for chunk in chunks)
    cooler.create_cooler(
        str(output_cooler_uri),
        bins[["chrom", "start", "end"]],
        chunks,
        dtypes={"count": float},
        ordered=True,
        assembly=assembly,
    )


def to_text_compressed(chunks, output_file: pathlib.Path, nproc: int):
    cmd = get_compressor(9, max(1, nproc - 1))
    with open(output_file, "wb") as f:
        if cmd is None:
            # Let pandas deal with compression (slow)
            return to_text(chunks, f)

        # For some reason using shell=True is faster
        with sp.Popen(cmd, stdin=sp.PIPE, stderr=sp.PIPE, stdout=f, shell=True) as compressor:
            to_text(chunks, compressor.stdin)
            compressor.communicate()
            if (code := compressor.returncode) != 0:
                print(compressor.stderr, file=sys.stderr)
                raise RuntimeError(f"{cmd} terminated with code {code}")


def to_text(chunks, output_file):
    if isinstance(output_file, str) or isinstance(output_file, pathlib.Path):
        with open(output_file, "wb") as f:
            return to_text(chunks, f)

    for chunk in chunks:
        chunk.drop(columns=["count"]).to_csv(
            output_file,
            sep="\t",
            index=False,
            header=False,
            na_rep="nan",
        )


def cooler_dump_chunked(
    cool_uri: str,
    weight_dset: str,
    join: bool,
    weight_type: str,
    chunksize=1000000,
):
    # Adapted from: https://github.com/open2c/cooler/blob/master/cooler/cli/dump.py
    clr = cooler.Cooler(cool_uri)

    # Load all the bins
    bins = clr.bins()[:]
    n_bins = len(bins)
    if chunksize is None:
        chunksize = len(bins)

    if weight_dset not in bins.columns:
        print("Balancing weights not found", file=sys.stderr)
        sys.exit(1)

    h5 = clr.open("r")
    reader = CSRReader(h5["pixels"], h5["indexes/bin1_offset"][:])
    field = "count"

    # Dump everything
    bbox = (0, n_bins, 0, n_bins)
    engine = DirectRangeQuery2D(reader, field, bbox, chunksize)

    chunks = (
        pd.DataFrame(
            dct,
            columns=["bin1_id", "bin2_id", field],
        )
        for dct in engine
    )

    if weight_type == "auto":
        divisive_weights = weights_are_divisive(weight_dset)
    else:
        divisive_weights = weight_type == "divisive"

    annotator = make_annotator(bins, weight_dset, join, divisive_weights)
    return map(annotator, chunks)


def main():
    args = vars(make_cli().parse_args())
    input_cooler = str(args["cooler"])
    output_file = args["output_name"]

    chunks = cooler_dump_chunked(
        input_cooler,
        args["norm_name"],
        output_file.suffix != ".cool",
        args["weight_type"],
    )

    if output_file.suffix == ".cool":
        bins = cooler.Cooler(input_cooler).bins()[:]
        to_cooler(chunks, bins, output_cooler_uri=args["output_name"])
    elif output_file.suffix == ".gz":
        to_text_compressed(chunks, output_file, args["nproc"])
    else:
        to_text(chunks, output_file)


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)

    for h in logging.getLogger("cooler").handlers:
        h.setFormatter(logging.Formatter(fmt))

    logging.getLogger("cooler.create").setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
