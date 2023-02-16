#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import itertools
import logging
import multiprocessing as mp
import os
import pathlib
import shlex
import shutil
import stat
import subprocess as sp
import sys
import tempfile
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

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    def find_executable(arg):
        if (cmd := shutil.which(arg)) is not None:
            return pathlib.Path(cmd)

        raise FileNotFoundError(f'Unable to find executable "{arg}"')

    supported_norms = {
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
        "ICE",
        "INTER_ICE",
        "GW_ICE",
    }

    cli.add_argument(
        "multires-cooler-in",
        type=pathlib.Path,
        help="Path to a .mcool file to use as input.",
    )

    cli.add_argument(
        "multires-cooler-out",
        type=pathlib.Path,
        help="Path to a .mcool file to use as output.",
    )

    cli.add_argument("--hic-tools-jar", type=existing_file, help="Path to the hic_tools jar file.")
    cli.add_argument("--juicer-tools-jar", type=existing_file, required=True, help="Path to the juicer_tools jar file.")
    cli.add_argument(
        "--hic2cool-ng", type=find_executable, default="hic2cool-ng", help="Path to hic2cool-ng executable."
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
        "--normalization-methods",
        type=str,
        nargs="+",
        choices=supported_norms,
        default=[norm for norm in supported_norms],
        help="One or more normalization methods.",
    )
    cli.add_argument(
        "--default-normalization-method",
        type=str,
        choices=supported_norms,
        default="ICE",
        help="Normalization method to store under weight/.",
    )
    cli.add_argument(
        "--tmpdir",
        type=pathlib.Path,
        default=pathlib.Path(tempfile.gettempdir()),
        help="Path to a folder to use for temporary file(s).",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    cli.add_argument(
        "-Xms",
        type=str,
        default="1024m",
    )
    cli.add_argument(
        "--compression-level",
        type=int,
        choices=range(1, 10),
        default=9,
        help="Gzip compression level used to compress temporary files.",
    )
    cli.add_argument(
        "-Xmx",
        type=str,
        default="32g",
    )
    cli.add_argument("--force", action="store_true", default=False, help="Overwrite existing files (if any).")

    return cli


def get_compressor(compression_level: int, nproc: int) -> Union[list, None]:
    if cmd := shutil.which("pigz"):
        return [cmd, f"-{compression_level}", "-p", str(min(4, nproc))]
    if cmd := shutil.which("gzip"):
        return [cmd, f"-{compression_level}"]
    return None


def read_blacklisted_regions(path_to_bed3: Union[None, List[pathlib.Path]]) -> Union[None, pd.DataFrame]:
    if path_to_bed3 is None or len(path_to_bed3) == 0:
        return pd.DataFrame({"chrom": [], "start": [], "end": []})

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = pd.concat([bf.read_table(p, schema="bed3") for p in path_to_bed3])

    df = bf.merge(df)[bf.SCHEMAS["bed3"]]
    logging.info("Imported %d regions for blacklisting", len(df))
    return df


def map_1d_regions_to_bin_ids(regions: pd.DataFrame, bins: pd.DataFrame) -> np.ndarray:
    if len(regions) == 0:
        return np.empty(0, dtype=int)
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


def fetch_pixels(cf: cooler.Cooler, chrom1: str, chrom2: str) -> pd.DataFrame:
    logging.debug("[%d] Processing pixels for %s:%s...", cf.binsize, chrom1, chrom2)
    pixels = pd.DataFrame(cf.matrix(balance=False, as_pixels=True, join=True).fetch(chrom1, chrom2))

    pixels["end1"] = 0
    pixels["end2"] = 1
    pixels.rename(columns={"end1": "str1", "end2": "str2"}, inplace=True)

    pixels["frag1"] = pixels["str1"]
    pixels["frag2"] = pixels["str2"]

    return pixels[["str1", "chrom1", "start1", "frag1", "str2", "chrom2", "start2", "frag2", "count"]]


def dump_pixels_helper_(uri: str, dest):
    num_pixels = 0
    cf = cooler.Cooler(uri)
    logging.info("[%d] Writing pixels to temporary file...", cf.binsize)

    queries = []
    chroms = cf.chromnames
    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i:]:
            queries.append([chrom1, chrom2])

    for chrom1, chrom2 in queries:
        pixels = fetch_pixels(cf, chrom1, chrom2)
        pixels.to_csv(dest, sep="\t", index=False, header=False)
        num_pixels += len(pixels)

    if isinstance(dest, pathlib.Path) or isinstance(dest, str):
        logging.info('[%s] Written a total of %s pixels to file "%s".', cf.binsize, num_pixels, dest)
    else:
        logging.info("[%d] Written a total of %s pixels.", cf.binsize, num_pixels)


def dump_pixels(uri: str, dest: pathlib.Path, compression_level: int, nproc: int) -> pathlib.Path:
    cmd = get_compressor(compression_level, max(1, nproc - 1))
    assert cmd is not None

    with open(dest, "wb") as f:  # For some reason using shell=True is faster
        with sp.Popen(cmd, stdin=sp.PIPE, stderr=sp.PIPE, stdout=f, shell=True) as compressor:
            dump_pixels_helper_(uri, compressor.stdin)
            compressor.communicate()
            if (code := compressor.returncode) != 0:
                print(compressor.stderr, file=sys.stderr)
                raise RuntimeError(f"{cmd} terminated with code {code}")

    return dest


def run_juicer_tools_pre(
    path_to_pixels: pathlib.Path,
    path_to_chrom_sizes: pathlib.Path,
    resolution: int,
    path_to_hic: pathlib.Path,
    path_to_jar: pathlib.Path,
    tmpdir: pathlib.Path,
    xms: str,
    xmx: str,
):
    logging.info("[%d] Converting matrix to .hic...", resolution)
    cmd = [
        shutil.which("java"),
        f"-Xms{xms}",
        f"-Xmx{xmx}",
        "-jar",
        str(path_to_jar),
        "pre",
        "-r",
        str(resolution),
        "-j",
        "2",
        "--threads",
        "2",
        "-t",
        str(tmpdir),
        "-n",
    ]

    cmd.extend([path_to_pixels, path_to_hic, path_to_chrom_sizes])

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.debug("[%d] Running %s...", resolution, cmd)
    sp.check_output(cmd, stderr=sp.STDOUT)


def cooler_to_hic(
    cooler_uri: str,
    hictools_jar: pathlib.Path,
    tmpdir: pathlib.Path,
    compression_lvl: int,
    nproc: int,
    xms: str,
    xmx: str,
) -> pathlib.Path:
    cf = cooler.Cooler(cooler_uri)
    with tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".chrom.sizes") as chromsizes, tempfile.NamedTemporaryFile(
        dir=tmpdir, suffix=".txt.gz"
    ) as pixels:
        cf.chromsizes.to_csv(chromsizes.name, sep="\t", index=True, header=False)
        dump_pixels(cooler_uri, pathlib.Path(pixels.name), compression_lvl, nproc)

        hic_file = tmpdir / f"{pathlib.Path(cf.filename).name}.{cf.binsize}.hic"

        run_juicer_tools_pre(
            pathlib.Path(pixels.name),
            pathlib.Path(chromsizes.name),
            cf.binsize,
            hic_file,
            hictools_jar,
            tmpdir,
            xms,
            xmx,
        )

        return hic_file


def run_juicer_tools_add_norm(
    path_to_hic: pathlib.Path,
    resolution: int,
    norm: str,
    path_to_jar: pathlib.Path,
    nproc: int,
    xms: str,
    xmx: str,
):

    if norm == "NONE":
        return

    logging.info("[%d] Normalizing using %s...", resolution, norm)

    cmd = [
        shutil.which("java"),
        f"-Xms{xms}",
        f"-Xmx{xmx}",
        "-jar",
        str(path_to_jar),
        "addNorm",
        "-j",
        str(nproc),
        "-F",
        "-w",
        str(resolution),
        "-r",
        str(resolution),
        "-k",
        norm,
        str(path_to_hic),
    ]

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.debug("[%s] Running %s...", resolution, cmd)
    sp.check_output(cmd, stderr=sp.STDOUT)


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
    cis_only = strategy == "ICE"
    trans_only = strategy == "INTER_ICE"

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


def write_weights(uri: str, strategy: str, bias: np.ndarray, stats: Dict):
    resolution = cooler.Cooler(uri).binsize
    logging.info("[%d] Writing weights to %s...", resolution, uri)
    fname = cooler.Cooler(uri).filename
    root = cooler.Cooler(uri).root.lstrip("/")

    path = f"{root}/bins"
    with h5py.File(fname, "r+") as cf:
        if strategy in cf[path]:
            del cf[path][strategy]

        h5opts = {"compression": "gzip", "compression_opts": 6}
        cf[path].create_dataset(strategy, data=bias, **h5opts)
        cf[path][strategy].attrs.update(stats)


def _juicer_tools_balance_helper(
    input_uri: str,
    output_uri: str,
    blacklist: pd.DataFrame,
    normalization_methods: List[str],
    hictools_jar: pathlib.Path,
    juicertools_jar: pathlib.Path,
    hic2cool_ng: pathlib.Path,
    tmpdir: pathlib.Path,
    compression_lvl: int,
    nproc: int,
    xms: str,
    xmx: str,
    lock,
):

    hic = None
    try:
        hic = cooler_to_hic(
            input_uri,
            hictools_jar,
            tmpdir,
            compression_lvl,
            nproc,
            xms,
            xmx,
        )

        resolution = cooler.Cooler(input_uri).binsize
        for norm in normalization_methods:
            run_juicer_tools_add_norm(hic, resolution, norm, juicertools_jar, nproc, xms, xmx)

            blacklist1d = map_1d_regions_to_bin_ids(blacklist, cooler.Cooler(input_uri).bins()[:])

            try:
                lock.acquire()
                run_hic2coolng_extract_norms(hic, output_uri, norm, hic2cool_ng)

                if len(blacklist1d) != 0:
                    path = cooler.Cooler(output_uri).filename
                    suffix = cooler.Cooler(output_uri).root
                    with h5py.File(path, "r+") as f:
                        weights = f[f"{suffix}/bins"][norm][:]
                        weights[blacklist1d] = np.nan

                        dset = f[f"{suffix}/bins/{norm}"]
                        dset[...] = weights

            finally:
                lock.release()

    finally:
        if hic is not None:
            hic.unlink()


def juicer_tools_balance(
    input_uri: str,
    output_uri: str,
    blacklist: pd.DataFrame,
    normalization_methods: List[str],
    hictools_jar: pathlib.Path,
    juicertools_jar: pathlib.Path,
    hic2cool_ng: pathlib.Path,
    tmpdir: pathlib.Path,
    compression_lvl: int,
    nproc: int,
    xms: str,
    xmx: str,
    pool,
    lock,
):
    input_uris = generate_multires_uri_list(input_uri)
    output_uris = generate_multires_uri_list(output_uri)

    pool.starmap(
        _juicer_tools_balance_helper,
        zip(
            input_uris,
            output_uris,
            itertools.repeat(blacklist),
            itertools.repeat(normalization_methods),
            itertools.repeat(hictools_jar),
            itertools.repeat(juicertools_jar),
            itertools.repeat(hic2cool_ng),
            itertools.repeat(tmpdir),
            itertools.repeat(compression_lvl),
            itertools.repeat(nproc),
            itertools.repeat(xms),
            itertools.repeat(xmx),
            itertools.repeat(lock),
        ),
        chunksize=1,
    )


def split_norm_strategies(strats: List[str]) -> Tuple[List[str], List[str]]:
    ice = []
    other = []

    for strat in strats:
        if strat.endswith("ICE"):
            ice.append(strat)
        else:
            other.append(strat)

    return ice, other


def run_hic2coolng_extract_norms(hic: pathlib.Path, uri: str, norm: str, hic2cool_ng: pathlib.Path):

    logging.info("[%d] Copying %s weights...", cooler.Cooler(uri).binsize, norm)
    cmd = [
        str(hic2cool_ng),
        "extract-norms",
        str(hic),
        uri,
        "--normalization-methods",
        norm,
        "--fail-if-norm-not-found",
    ]

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.debug("Running %s...", cmd)
    sp.check_output(cmd, stderr=sp.STDOUT)


def copy_weights(uri: str, src_weight_dset: str, dest_weight_dset: str):
    resolution = cooler.Cooler(uri).binsize
    logging.info("[%d] Copying %s to %s...", resolution, src_weight_dset, dest_weight_dset)
    fname = cooler.Cooler(uri).filename
    root = cooler.Cooler(uri).root.lstrip("/")

    path = f"{root}/bins"
    with h5py.File(fname, "r+") as cf:
        if dest_weight_dset in cf[path]:
            del cf[path][dest_weight_dset]

        src = f"{path}/{src_weight_dset}"
        cf.copy(cf[src], cf[path], dest_weight_dset)


def write_default_weight_dataset(path: pathlib.Path, default_norm: str):
    for uri in generate_multires_uri_list(str(path)):
        copy_weights(uri, default_norm, "weight")


def main():
    args = vars(make_cli().parse_args())

    default_norm = args["default_normalization_method"]
    norm_methods = args["normalization_methods"]
    if default_norm not in norm_methods:
        norm_methods = ", ".join(norm_methods)
        raise RuntimeError(
            f"{default_norm} not in {norm_methods}.\n"
            "Please specify a different default normalization method with --default-normalization-method."
        )

    input_cooler = str(args["multires-cooler-in"])
    output_cooler = str(args["multires-cooler-out"])

    blacklist = read_blacklisted_regions(args.get("blacklist"))

    juicer_tools_jar = args["juicer_tools_jar"]
    if (hic_tools_jar := args.get("hic_tools_jar")) is None:
        hic_tools_jar = juicer_tools_jar

    ice_norms, other_norms = split_norm_strategies(args["normalization_methods"])
    with mp.Pool(args["nproc"]) as pool, mp.Manager() as manager:
        with tempfile.TemporaryDirectory(dir=args["tmpdir"]) as tmpdir:
            if input_cooler == output_cooler:
                tmp_cooler = pathlib.Path(tmpdir) / input_cooler
                shutil.copy(input_cooler, tmp_cooler)
                input_cooler = tmp_cooler
            else:
                shutil.copy(input_cooler, output_cooler)
                os.chmod(output_cooler, os.stat(output_cooler).st_mode | stat.S_IWUSR)

            check_rw_permissions_on_cooler(output_cooler)
            if len(other_norms) != 0:
                juicer_tools_balance(
                    input_cooler,
                    output_cooler,
                    blacklist=blacklist,
                    normalization_methods=other_norms,
                    hictools_jar=hic_tools_jar,
                    juicertools_jar=juicer_tools_jar,
                    hic2cool_ng=args["hic2cool_ng"],
                    tmpdir=pathlib.Path(tmpdir),
                    compression_lvl=args["compression_level"],
                    nproc=args["nproc"],
                    xms=args["Xms"],
                    xmx=args["Xmx"],
                    pool=pool,
                    # Protect write access to output cooler file
                    lock=manager.Lock(),
                )

        if len(ice_norms) != 0:
            for uri, strategy in generate_tasks(output_cooler, ice_norms):
                bias, stats = run_cooler_balance(
                    uri, strategy, blacklist, args["nproc"], pool, mad_max=args["mad_max_threshold"]
                )
                write_weights(uri, strategy, bias, stats)

        write_default_weight_dataset(pathlib.Path(output_cooler), default_norm)


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
