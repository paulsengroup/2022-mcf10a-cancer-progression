#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import hashlib
import pathlib
import random
import sys
from typing import List


def make_cli() -> argparse.ArgumentParser:
    def positive_int(s) -> float:
        if (x := int(s)) > 0:
            return x

        raise RuntimeError("Not a positive int")

    cli = argparse.ArgumentParser(
        "Use one or more file as source of entropy to generate a sequence of seeds. Seeds are printed to stdout."
    )
    cli.add_argument(
        "files",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more files to use as source of entropy",
    )

    cli.add_argument(
        "--number-of-seeds",
        type=positive_int,
        default=50,
        help="Number of seeds to generate.",
    )
    cli.add_argument("--sep", type=str, default="\n", help="Seed separator.")
    cli.add_argument(
        "--lower-bound",
        type=int,
        default=0,
        help="Lower bound for the uniform distribution from which seeds are sampled.",
    )
    cli.add_argument(
        "--upper-bound",
        type=int,
        default=(2**32) - 1,
        help="Upper bound for the uniform distribution from which seeds are sampled.",
    )

    return cli


def hash_files(files: List[pathlib.Path], chunk_size: int = 64 * 1024 * 1024) -> str:
    hasher = hashlib.sha512()

    for f in set(files):
        with open(f, "rb") as fp:
            while data := fp.read(chunk_size):
                hasher.update(data)

    return hasher.hexdigest()


def main():
    args = vars(make_cli().parse_args())

    digest = hash_files(args["files"])

    random.seed(digest, version=2)

    lb = args["lower_bound"]
    ub = args["upper_bound"]

    num_seeds = args["number_of_seeds"]
    num_seeds_str_length = len(str(num_seeds))

    seeds = (f"{i:0{num_seeds_str_length}d}\t{random.randint(lb, ub)}" for i in range(num_seeds))
    print(args["sep"].join(seeds))


if __name__ == "__main__":
    main()
