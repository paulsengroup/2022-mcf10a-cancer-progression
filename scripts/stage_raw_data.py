#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
from collections import namedtuple


def parse_file_name(path, pattern):
    res = pattern.search(path)
    res = res.groups()
    assert len(res) == 3

    MatchT = namedtuple("MatchT", ["id", "name", "read_pair"])
    return MatchT(int(res[0]), res[1], int(res[2]))


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("--input-folder", type=str, required=True)
    cli.add_argument("--output-folder", type=str, required=True)
    cli.add_argument(
        "--dry-run", action="store_const", default=True, const=True, dest="dry_run"
    )
    cli.add_argument("--no-dry-run", action="store_const", const=False, dest="dry_run")
    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()
    sample_names = tuple(
        ["10Anew", "10A", "10_alpha", "T1", "C1", "TGFB", "shHPA_Z", "shZ"]
    )
    p = re.compile(
        r"^.*/HiC_(\d+)_(" + "|".join(sample_names) + r").*_S\d+_R([12]).*.fastq.gz$"
    )

    print(f"mkdir -p {args.output_folder}", file=sys.stderr)
    if not args.dry_run:
        os.makedirs(args.output_folder, exist_ok=True)

    for inpath in glob.glob(f"{args.input_folder}/*.fastq.gz"):
        toks = parse_file_name(inpath, p)
        outpath = f"{args.output_folder}/HiC_{toks.id:03d}_{toks.name}_R{toks.read_pair}.fastq.gz"
        print(f'ln -s "{inpath}" "{outpath}"', file=sys.stderr)
        if not args.dry_run:
            os.symlink(inpath, outpath)
