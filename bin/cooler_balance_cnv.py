#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import os
from typing import Dict, List, Tuple

import hictkpy as htk
import iced
import numpy as np
import numpy.typing as npt
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "input-cooler",
        type=str,
        help="Path to a .cool file to be processed (URI syntax is supported).",
    )
    cli.add_argument(
        "cnv-bedgraph",
        type=str,
        help="Path to a Bedgraph file with CNV information.",
    )
    cli.add_argument(
        "output-cooler",
        type=str,
        help="Path to the output .cool file.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing file(s).",
    )
    cli.add_argument(
        "--chromosomes",
        type=str,
        help="Comma-separated list of chromosomes to be processed.",
    )
    cli.add_argument(
        "--method",
        type=str,
        choices={"CAIC", "LOIC"},
        required=True,
        help="Balancing method.",
    )

    return cli


def map_old_bins_to_new_bins(bins: pd.DataFrame, chroms: List[str]) -> Dict[int, int]:
    bins = (
        bins[bins["chrom"].isin(chroms)]
        .reset_index()
        .rename(columns={"index": "old_index"})
        .reset_index()
        .rename(columns={"index": "new_index"})
        .drop(columns={"chrom", "start", "end"})
    )

    bin_map = {}
    for new_index, old_index in bins.itertuples(index=False):
        bin_map[old_index] = new_index

    return bin_map


def map_new_bins_to_old_bins(bins: pd.DataFrame, chroms: List[str]) -> Dict[int, int]:
    bins = (
        bins[bins["chrom"].isin(chroms)]
        .reset_index()
        .rename(columns={"index": "old_index"})
        .reset_index()
        .rename(columns={"index": "new_index"})
        .drop(columns={"chrom", "start", "end"})
    )

    bin_map = {}
    for new_index, old_index in bins.itertuples(index=False):
        bin_map[new_index] = old_index

    return bin_map


def build_matrix(clr: htk.File, chroms: List[str]) -> Tuple[npt.NDArray, Dict[int, int]]:
    bin_map = map_old_bins_to_new_bins(clr.bins(), chroms)

    m = np.zeros([len(bin_map)] * 2)

    for i1, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i1:]:
            for pxl in clr.fetch(chrom1, chrom2):
                bin1 = bin_map[pxl.bin1_id]
                bin2 = bin_map[pxl.bin2_id]
                m[bin1, bin2] = pxl.count
    return m, map_new_bins_to_old_bins(clr.bins(), chroms)


def numpy2d_to_cooler(
    output_path: str, m: npt.NDArray, bin_mappings: Dict[int, int], chromosomes: Dict[str, int], resolution: int
):
    data = []
    idx = np.where(m != 0)
    for i1, i2 in zip(idx[0], idx[1]):
        data.append([bin_mappings[i1], bin_mappings[i2], m[i1, i2]])

    df = pd.DataFrame(data, columns=["bin1_id", "bin2_id", "count"])

    assert (df["count"] == 0).sum() == 0

    writer = htk.cooler.FileWriter(output_path, chromosomes, resolution)
    writer.add_pixels(df)
    writer.finalize()


def main():
    args = vars(make_cli().parse_args())

    output_path = args["output-cooler"]
    if not args["force"] and os.path.exists(output_path):
        raise RuntimeError(f'Refusing to overwrite output file "{output_path}"')
    else:
        try:
            os.remove(args["output-cooler"])
        except FileNotFoundError:
            pass

    cnv = pd.read_table(args["cnv-bedgraph"], names=["chrom", "start", "end", "cnv"], usecols=list(range(4)))
    clr = htk.File(args["input-cooler"], 0)

    if args["chromosomes"] is None:
        m = np.triu(clr.fetch().to_numpy()).astype(float)
        cnv = cnv["cnv"].to_numpy()
        bin_mappings = {i: i for i in range(len(clr.bins()))}
    else:
        chroms = args["chromosomes"].split(",")
        m, bin_mappings = build_matrix(clr, chroms)
        cnv = cnv.loc[cnv["chrom"].isin(chroms), "cnv"].to_numpy()

    assert len(cnv) == len(m)
    m = iced.normalization.ICE_normalization(m, counts_profile=cnv)

    if args["method"] == "CAIC":
        bins = clr.bins()
        if args["chromosomes"] is not None:
            chroms = args["chromosomes"].split(",")
            bins = bins[bins["chrom"].isin(chroms)]

        num_bins = bins.groupby("chrom", observed=True).size().to_numpy()
        block_biases = iced.normalization.estimate_block_biases(m, num_bins, cnv)
        m /= block_biases

    numpy2d_to_cooler(output_path, m, bin_mappings, clr.chromosomes(), clr.resolution())


if __name__ == "__main__":
    main()
