# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.0.0 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG COOLTOOLS_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

# ucsc-bedgraphtobigwig is required by certain cooltools subcommands
RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        python=3.10 \
        r=4.2 \
        bioframe \
        cooler \
        matplotlib \
        numpy \
        pandas \
        pyBigWig \
        r-alluvial \
        r-data.table \
        r-optparse \
        r-stringr \
        r-svglite \
        seaborn \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"

ENTRYPOINT ["/bin/bash"]
WORKDIR /data


LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-plotting}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

