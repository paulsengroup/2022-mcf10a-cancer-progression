# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.4.3 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        python=3.10 \
        'bioframe>=0.4.1' \
        'cooler>=0.9.2' \
        'hictkpy==0.0.5' \
        matplotlib \
        matplotlib-venn \
        natsort \
        numpy \
        'pandas>=2' \
        procps-ng \
        pyBigWig \
        seaborn \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data


LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-pyutils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
