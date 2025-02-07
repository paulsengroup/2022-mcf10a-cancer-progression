#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


# The point of this script is to compare the original FASTQ files with the ones downloaded from SRA
# to make sure they are identical (except for the read identifiers, which are modified by SRA
# upon submission...)


import argparse
import gzip
import hashlib
import io
import itertools
import logging
import multiprocessing as mp
import pathlib
import random
import time

import pandas as pd
import requests


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to the download.tsv file with the list of SRA FASTQ files to be compared.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        default=1,
    )

    return cli


def download_and_compare(row: pd.Series):
    url = row[1]["url"]
    path = row[1]["dest"]

    logging.info("begin processing %s...", url)
    status_code = None
    response = None
    for i in range(5):
        sleep_time = random.uniform(0, 180)
        logging.info("sleeping for %.2f seconds...", sleep_time)
        time.sleep(sleep_time)

        logging.info("downloading %s...", url)

        response = requests.get(url, stream=True)

        if response.status_code == 200:
            break

        logging.warning("failed to download %s (attempt %d/5), retrying...", url, i + 1)
        status_code = response.status_code
        response = None

    if response is None:
        raise RuntimeError(f'Failed to download file "{url}": status code: {status_code}')

    hasher1 = hashlib.sha256()
    hasher2 = hashlib.sha256()

    with gzip.GzipFile(fileobj=response.raw, mode="rb") as gz1, gzip.GzipFile(path, mode="rb") as gz2:
        with io.TextIOWrapper(gz1, encoding="utf-8") as f1, io.TextIOWrapper(gz2, encoding="utf-8") as f2:
            logging.info("beginning to parse %s and comparing its content with %s...", url, path)
            for i, (line1, line2) in enumerate(itertools.zip_longest(f1, f2)):
                if line1 is None or line2 is None:
                    raise RuntimeError(f"{url}: file lengths differ!")
                if line1.startswith("@"):
                    assert line2.startswith("@")
                    continue
                hasher1.update(line1.encode("utf-8"))
                hasher2.update(line2.encode("utf-8"))

                if (i + 1) % 100_000_000 == 0:
                    digest1 = hasher1.hexdigest()
                    digest2 = hasher2.hexdigest()
                    if digest1 == digest2:
                        logging.info("%s[:%d]: OK!", url, i + 1)
                    else:
                        logging.error("%s[:%d]: FAILURE! Expected %s, found %s", url, i + 1, digest2, digest1)
                        raise RuntimeError(f'URL "{url}": expected {digest1}, found {digest2}')

    digest1 = hasher1.hexdigest()
    digest2 = hasher2.hexdigest()

    if digest1 == digest2:
        logging.info("%s: File successfully validated!", url)
    else:
        logging.error("%s: FAILURE! Expected %s, found %s", url, digest2, digest1)
        raise RuntimeError(f'URL "{url}": expected {digest1}, found {digest2}')


def setup_logger(level: str):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=fmt)
    logging.getLogger().setLevel(level)


def main():
    setup_logger("INFO")
    args = vars(make_cli().parse_args())

    df = pd.read_table(args["tsv"])
    df["url"] = df["url"].str.replace("^ftp://", "https://", regex=True)
    df = df[df["url"].str.startswith("https://ftp.sra.ebi.ac.uk")]

    with mp.Pool(args["nproc"]) as pool:
        pool.map(download_and_compare, df.iterrows())


if __name__ == "__main__":
    main()
