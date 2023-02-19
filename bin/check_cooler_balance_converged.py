#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
from typing import List, Union

import cooler
import h5py
import numpy as np


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
        type=str,
        help="Path to a .cool or .mcool file to validate.",
    )
    cli.add_argument("--normalization-dataset", type=str, default="weight")

    return cli


def read_convergence_attribure(uri: str, dset_name: str) -> Union[None, bool, np.ndarray]:
    path, suffix = parse_cooler_uri(uri)
    with h5py.File(path) as h5:
        attrs = h5[f"{suffix}/bins/{dset_name}"].attrs
        if (converged := attrs.get("converged")) is not None:
            return converged

    return None


def generate_uri_list(uri: str) -> List[str]:
    if cooler.fileops.is_cooler(uri):
        return [uri]

    uris = []
    for suffix in cooler.fileops.list_coolers(uri):
        uris.append(f"{uri}::{suffix}")

    return uris


def main():
    args = vars(make_cli().parse_args())

    all_converged = True
    for uri in generate_uri_list(args["cooler"]):
        converged = read_convergence_attribure(uri, args["normalization_dataset"])
        if converged is None:
            logging.warning("%s does not have the convergence attribute", uri)
        elif isinstance(converged, bool):
            lvl = logging.INFO if converged else logging.WARNING
            logging.log(lvl, "%s: converged: %s", uri, converged)
            all_converged &= converged
        else:
            converged_ = np.all(converged)
            lvl = logging.INFO if converged_ else logging.WARNING
            logging.log(lvl, "%s: converged: %s", uri, converged_)
            all_converged &= converged_
            if not converged_:
                chroms = np.array(cooler.Cooler(uri).chromnames)
                for chrom in chroms[~converged]:
                    logging.info("%s: %s did not converge", uri, chrom)

    print(f"converged: {all_converged}")

    return not all_converged


def setup_logger(level=logging.WARNING):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
