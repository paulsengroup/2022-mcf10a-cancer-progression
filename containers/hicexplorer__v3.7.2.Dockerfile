# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:0.27.0 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG HICEXPLORER_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

# ucsc-bedgraphtobigwig is required by certain cooltools subcommands
RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        "hicexplorer=$HICEXPLORER_VER" \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /data


RUN hicConvertFormat --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-david-hic'
LABEL org.opencontainers.image.documentation='https://github.com/2022-david-hic'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-david-hic'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hicexplorer}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

