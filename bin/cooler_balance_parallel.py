#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import itertools
import logging
import multiprocessing as mp
import pathlib
import warnings
from typing import Dict, List, Tuple, Union

import bioframe as bf
import cooler
import h5py
import numpy as np
import pandas as pd


def parse_cooler_uri(uri):
    cf = cooler.Cooler(uri)
    return cf.filename, cf.root.lstrip("/")


def printable_chrom(chrom):
    if chrom is None:
        return "all"
    return str(chrom)


def make_cli():
    cli = argparse.ArgumentParser()

    def positive_float(s) -> float:
        if (x := float(s)) > 0:
            return x

        raise RuntimeError("Not a positive float")

    cli.add_argument(
        "cooler",
        type=pathlib.Path,
        help="Path to a .cool or .mcool file (URI syntax supported).\n"
        "When passing an .mcool, all resolutions will be balanced.",
    )
    cli.add_argument(
        "--blacklist",
        type=pathlib.Path,
        nargs="+",
        help="Path to one or more a BED3+ files with the list of regions to mask out during blancing.",
    )
    cli.add_argument(
        "--mad-max-threshold",
        type=positive_float,
        default=10.0,
        help="Threshold for the MAD-max filter (see cooler's docs for more details).",
    )
    cli.add_argument(
        "--strategy",
        type=str,
        nargs="+",
        choices={"genomewide", "cis-only", "trans-only"},
        default=["genomewide", "cis-only", "trans-only"],
        help="One or more normalization strategy.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )

    return cli


def read_blacklisted_regions(path_to_bed3: Union[None, List[pathlib.Path]]) -> Union[None, pd.DataFrame]:
    if path_to_bed3 is None or len(path_to_bed3) == 0:
        return None

    with warnings.catch_warnings():
        # Suppress UserWarning regarding unstable API
        warnings.simplefilter("ignore")
        df = pd.concat([bf.read_table(p, schema="bed3") for p in path_to_bed3])

    df = bf.merge(df)[bf.SCHEMAS["bed3"]]
    logging.info("Imported %d regions for blacklisting", len(df))
    return df


def map_1d_regions_to_bin_ids(regions: pd.DataFrame, bins: pd.DataFrame) -> np.ndarray:
    df = bf.overlap(regions, bins[bf.SCHEMAS["bed3"]], return_index=True, suffixes=("_1", "_2")).dropna()
    return df["index_2"].to_numpy().astype(int)


def check_rw_permissions_on_cooler(uri: str):
    if cooler.fileops.is_multires_file(uri):
        suffix = cooler.fileops.list_coolers(uri)[0]
        uri = f"{uri}::{suffix}"

    path, suffix = parse_cooler_uri(uri)
    with h5py.File(path, "r+") as f:
        pass


def generate_multires_uri_list(uri: str) -> List[str]:
    uris = []
    for suffix in cooler.fileops.list_coolers(uri):
        uris.append(f"{uri}::{suffix}")

    return uris


def generate_tasks(uri: str, strategies: List[str]) -> List[Tuple[str, str]]:
    if cooler.fileops.is_multires_file(uri):
        uris = generate_multires_uri_list(uri)
    else:
        uris = [uri]

    return [tuple(x) for x in itertools.product(uris, strategies)]


# Set param default values to the same values used by cooler balance from cooler v0.9.1.
# This is done to avoid issues when older/newer versions of cooler are used
def run_cooler_balance(
    uri: str,
    strategy: str,
    blacklist: Union[None, pd.DataFrame],
    num_chunks: int,
    pool,
    ignore_diags=2,
    mad_max=5,
    min_nnz=10,
    min_count=0,
    rescale_marginals=True,
    tol=1e-05,
    max_iters=2000,
) -> Tuple[np.ndarray, Dict]:
    cis_only = strategy == "cis-only"
    trans_only = strategy == "trans-only"

    cf = cooler.Cooler(uri)
    bin_size = cf.binsize
    num_bins = cf.info["nbins"]

    if blacklist is None:
        blacklist1d = None
    else:
        blacklist1d = map_1d_regions_to_bin_ids(blacklist, cf.bins()[:])

        logging.info(
            "[%d - %s] blacklisting %d bins (%.2g%%)",
            bin_size,
            strategy,
            len(blacklist1d),
            100 * len(blacklist1d) / num_bins,
        )

    if cis_only:
        chunk_size = 5000000
        logging.info("[%d - %s] balancing using chunks of size %d...", bin_size, strategy, chunk_size)
    else:
        chunk_size = int(cf.info["nnz"] / num_chunks)
        logging.info("[%d - %s] balancing using %d chunks...", bin_size, strategy, num_chunks)

    return cooler.balance_cooler(
        cf,
        chunksize=chunk_size,
        cis_only=cis_only,
        trans_only=trans_only,
        blacklist=blacklist1d,
        map=pool.map,
        ignore_diags=ignore_diags,
        mad_max=mad_max,
        min_nnz=min_nnz,
        min_count=min_count,
        rescale_marginals=rescale_marginals,
        tol=tol,
        max_iters=max_iters,
        store=False,
    )


def write_weights(uri: str, strategy: str, bias: np.ndarray, stats: Dict, store_name: str = "weight"):
    logging.info("Writing weights to %s...", uri)
    fname = cooler.Cooler(uri).filename
    root = cooler.Cooler(uri).root.lstrip("/")
    if strategy != "genomewide":
        store_name += f"-{strategy}"

    path = f"{root}/bins"
    with h5py.File(fname, "r+") as cf:
        if store_name in cf[path]:
            del cf[path][store_name]

        h5opts = {"compression": "gzip", "compression_opts": 6}
        cf[path].create_dataset(store_name, data=bias, **h5opts)
        cf[path][store_name].attrs.update(stats)


def main():
    args = vars(make_cli().parse_args())
    input_cooler = str(args["cooler"])

    check_rw_permissions_on_cooler(input_cooler)

    blacklist = read_blacklisted_regions(args.get("blacklist"))

    with mp.Pool(args["nproc"]) as pool:
        for uri, strategy in generate_tasks(input_cooler, args["strategy"]):
            bias, stats = run_cooler_balance(
                uri, strategy, blacklist, args["nproc"], pool, mad_max=args["mad_max_threshold"]
            )
            write_weights(uri, strategy, bias, stats)


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)

    for h in logging.getLogger("cooler").handlers:
        h.setFormatter(logging.Formatter(fmt))

    logging.getLogger("cooler.balance").setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
