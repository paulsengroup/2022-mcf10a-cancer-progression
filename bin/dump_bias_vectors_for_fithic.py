#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import cooler
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser(
        "Output the bias vector required by dcHiC/FitHiC.\n"
        "The script assumes weights have been computed using cooler, meaning that:\n"
        " - cooler file is expected to have a scale attribute\n"
        " - weights should be multiplicative"
    )

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
    return cli


def read_bins(clr: cooler.Cooler, norm_name: str) -> pd.DataFrame:
    return clr.bins()[["chrom", "start", "end", norm_name]][:].rename(columns={norm_name: "weight"})


def read_weight_scale_attribute(clr: cooler.Cooler, weight_dset_name: str) -> pd.DataFrame:
    scale = np.sqrt(clr.open("r")[f"bins/{weight_dset_name}"].attrs["scale"])
    if isinstance(scale, float):
        return pd.DataFrame(data={"chrom": clr.chromnames, "scale": [scale] * len(clr.chromnames)})

    return pd.DataFrame(data={"chrom": clr.chromnames, "scale": scale})


def std_scaler(values: pd.Series) -> pd.Series:
    mean = values.mean(skipna=True, numeric_only=True)
    std = values.std(skipna=True, numeric_only=True)

    return (values - mean) / std


def scale_weights(bins: pd.DataFrame) -> pd.Series:
    """
    Center weights around 1
    """
    return 1.0 + bins.groupby("chrom").transform(std_scaler)["weight"]


def process_bins(bins: pd.DataFrame, scale: pd.DataFrame) -> pd.DataFrame:
    bins = bins.merge(scale, on="chrom", how="left")

    bin_size = (bins["end"] - bins["start"]).max()

    # This often leads to position past the end of chromosomes when processing the last bin
    # However this is what FitHiC expects (as in, positions are expected to always be multple of bin_size / 2)
    bins["pos"] = bins["start"] + int(round(bin_size / 2))

    bins["weight"] *= bins["scale"]
    bins["weight"] = scale_weights(bins)
    bins["weight"] = np.nan_to_num(bins["weight"], nan=-1)

    return bins[["chrom", "pos", "weight"]]


def main():
    args = vars(make_cli().parse_args())

    clr = cooler.Cooler(str(args["cooler"]))
    bins = process_bins(
        read_bins(clr, args["norm_name"]),
        read_weight_scale_attribute(clr, args["norm_name"]),
    )
    bins.to_csv(sys.stdout, index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()
