#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import multiprocessing as mp
import pathlib
import sys
from typing import List, Union

import cooler
import hicrep
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    def existing_cooler(arg):
        return cooler.Cooler(arg).uri

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli.add_argument(
        "coolers",
        nargs="+",
        type=existing_cooler,
        help="Path to 2 or more Cooler files (URI syntax supported).",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma separated list of sample labels to use instead of labels inferred from file names.",
    )

    cli.add_argument(
        "--h",
        type=positive_int,
        required=True,
        help="Smooth the input contact matrices using a 2d mean filter with window size of 1 + 2 * value. "
        "This should be set according to the bin size.",
    )

    cli.add_argument(
        "--dBPMax",
        type=positive_int,
        default=int(5e6),
        help="Only consider contacts at most this number of bp away from the diagonal.",
    )

    cli.add_argument(
        "--bDownSample",
        action="store_true",
        default=False,
        help="Down sample the input with more contact counts to the the same number of counts as the other input "
        "with less contact counts.",
    )

    cli.add_argument(
        "--nproc",
        "-p",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=1,
        help="Maximum number of CPU cores to use for parallel processing.",
    )

    return cli


def validate_coolers(cooler_uris: List[str], labels: Union[List[str], None]) -> None:
    # TODO: validate chromosome names and sizes
    if len(cooler_uris) < 2:
        raise RuntimeError(f"Expected 2 or more Coolers, found {len(cooler_uris)}.")

    if labels is not None and len(labels) != len(cooler_uris):
        raise RuntimeError(
            f"Mismatch in the number of Coolers and labels, found {len(cooler_uris)} and {len(labels)} respectively."
        )

    bin_sizes = set()
    for uri in cooler_uris:
        c = cooler.Cooler(uri)
        bin_sizes.add(c.binsize)

    if len(bin_sizes) != 1:
        bin_sizes = ", ".join((str(bs) for bs in bin_sizes))
        raise RuntimeError(
            f"URIs point to Coolers with different resolutions. Found the following resolutions: {bin_sizes}"
        )


def run_hicrep(
    cooler_uri1: str,
    cooler_uri2: str,
    label1: Union[str, None],
    label2: Union[str, None],
    h: int,
    dbpmax: int,
    downsample: bool,
) -> pd.DataFrame:
    cooler1 = cooler.Cooler(cooler_uri1)
    cooler2 = cooler.Cooler(cooler_uri2)

    num_chroms = len(cooler1.chromsizes)
    if label1 is None:
        assert label2 is None
        labels1 = [pathlib.Path(cooler1.filename).name] * num_chroms
        labels2 = [pathlib.Path(cooler2.filename).name] * num_chroms
    else:
        labels1 = [label1] * num_chroms
        labels2 = [label2] * num_chroms

    if cooler1.uri == cooler2.uri:
        chroms = cooler1.chromnames
        scc = [1.0] * num_chroms

        return pd.DataFrame(
            {"cond1": labels1, "cond2": labels2, "chrom": chroms, "scc": scc}
        )

    assert (cooler1.chromsizes == cooler2.chromsizes).all()

    scc = hicrep.hicrepSCC(cooler1, cooler2, h, dbpmax, downsample)

    assert len(scc) == num_chroms

    return pd.DataFrame(
        {
            "cond1": labels1,
            "cond2": labels2,
            "chrom": cooler1.chromnames,
            "scc": scc,
        }
    )


def generate_hicrep_args(
    coolers: List[str],
    labels: Union[List[str], None],
    h: int,
    dbpmax: int,
    downsample: bool,
) -> list:
    args = []
    for i in range(len(coolers)):
        for j in range(i, len(coolers)):
            arg = [coolers[i], coolers[j]]
            if labels is None:
                arg.extend([None, None])
            else:
                arg.extend([labels[i], labels[j]])

            args.append(arg + [h, dbpmax, downsample])

    return args


def main():
    args = vars(make_cli().parse_args())

    coolers = args["coolers"]
    labels = args.get("labels")
    if labels is not None:
        labels = labels.split(",")
    validate_coolers(coolers, labels)

    hicrep_args = generate_hicrep_args(
        coolers, labels, args["h"], args["dBPMax"], args["bDownSample"]
    )

    with mp.Pool(args["nproc"]) as pool:
        dfs = pool.starmap(
            run_hicrep,
            hicrep_args,
            chunksize=1,
        )

    pd.concat(dfs).to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
