#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import hashlib
import logging
import multiprocessing as mp
import pathlib
import random
import sys
import time

import pandas as pd
import requests


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to the download.tsv file with the list of SRA FASTQ files to be checksummed.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        default=1,
    )

    return cli


def download_and_checksum(row: pd.Series) -> pd.Series:
    row = row[1].copy()
    url = row["url"]

    if not url.startswith("https://ftp.sra.ebi.ac.uk"):
        return row

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

    hasher = hashlib.sha256()

    while True:
        data = response.raw.read(64 << 20)
        if not data:
            break
        hasher.update(data)

    row["sha256"] = hasher.hexdigest()

    return row


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
        df = pd.DataFrame.from_records(pool.map(download_and_checksum, df.iterrows())).sort_values("dest")
        df["url"] = df["url"].str.replace("^https://", "ftp://", regex=True)
        df.to_csv(sys.stdout, index=False, sep="\t")


if __name__ == "__main__":
    main()
