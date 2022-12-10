# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.1.0 AS downloader
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_NO_CACHE_DIR=0

ARG CONTAINER_VERSION
ARG MAXHIC_VER=${CONTAINER_VERSION}

RUN micromamba install -y -c conda-forge git \
&& git clone https://github.com/bcb-sut/MaxHiC /tmp/maxhic \
&& cd /tmp/maxhic \
&& git checkout "${MAXHIC_VER}"

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

FROM mambaorg/micromamba:1.1.0 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG MAXHIC_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        cooler \
        "numpy>=1.4,<2" \
        "pandas>=0.24,<1" \
        "scipy>=1.1,<2" \
        "tensorflow>=1.3,<2" \
&& micromamba clean --all -y


COPY --from=downloader --chown=nobody:nogroup /tmp/maxhic/Capture /opt/maxhic/Capture
COPY --from=downloader --chown=nobody:nogroup /tmp/maxhic/General /opt/maxhic/General
COPY --from=downloader --chown=nobody:nogroup /tmp/maxhic/Main.py /opt/maxhic/

USER root
RUN mkdir /opt/maxhic/bin/ \
&& printf '%s\npython /opt/maxhic/Main.py "$@"\n' '#!/bin/sh' > /opt/maxhic/bin/maxhic \
&& chmod 755 /opt/maxhic/bin/maxhic \
&& chown -R nobody:nogroup /opt/maxhic/
USER mambauser

ENV PATH="/opt/conda/bin:/opt/maxhic/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/maxhic/bin/maxhic"]

WORKDIR /data

RUN maxhic --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-maxhic}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
