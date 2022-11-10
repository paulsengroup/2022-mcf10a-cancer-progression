# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.0.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DESEQ2_VERSION="1.38.*"

RUN micromamba install -y                            \
               -c conda-forge                        \
               -c bioconda                           \
               "bioconductor-deseq2=$DESEQ2_VERSION" \
               bioframe                              \
&& micromamba clean --all -y

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-david-hic'
LABEL org.opencontainers.image.documentation='https://github.com/2022-david-hic'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-david-hic'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-diff-expression-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
