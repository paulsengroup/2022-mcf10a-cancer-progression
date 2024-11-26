# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:2.0.4 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG CLUSTER_PROFILER="$CONTAINER_VERSION"
ARG ORG_HS_EG_DB_VER=3.17

RUN micromamba install -y                                       \
               -c conda-forge                                   \
               -c bioconda                                      \
               python                                           \
               r                                                \
               "bioconductor-clusterprofiler=$CLUSTER_PROFILER" \
               "bioconductor-org.hs.eg.db=$ORG_HS_EG_DB_VER"    \
               bioframe                                         \
               numpy                                            \
               pandas                                           \
               procps-ng                                        \
               rpy2                                             \
&& micromamba clean --all -y


ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-cluster-profiler}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
