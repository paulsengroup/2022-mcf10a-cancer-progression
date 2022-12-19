# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG HICREP_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
&& apt-get install -y -q --no-install-recommends \
                      cython3 \
                      gcc \
                      python3 \
                      python3-dev \
                      python3-pip \
&& pip install "hicrep==$HICREP_VER" \
               "matplotlib<3.6" \
               "numpy<1.24" \
&& apt-get remove -y -q cython3 \
                         gcc \
                         python3-dev \
                         python3-pip \
&& apt-get autoremove -y \
&& rm -rf /var/lib/apt/lists/*

CMD ["/usr/local/bin/hicrep"]
WORKDIR /data

RUN hicrep --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hicrep}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

