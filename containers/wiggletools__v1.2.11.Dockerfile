# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.5.10 AS base

ARG CONTAINER_VERSION
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y                         \
               -c conda-forge                     \
               -c bioconda                        \
               procps-ng                          \
               ucsc-wigtobigwig                   \
               "wiggletools=${CONTAINER_VERSION}" \
&& micromamba clean --all -y


WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

RUN wiggletools --help
RUN whereis wigToBigWig

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-wiggletools}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
