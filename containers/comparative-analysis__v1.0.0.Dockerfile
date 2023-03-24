# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.4.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y \
               -c conda-forge \
               -c bioconda \
               bioframe \
               matplotlib \
               numpy \
               pandas \
               pigz \
               pyBigWig \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

RUN python3 -c 'import bioframe, matplotlib, numpy, pandas, pyBigWig'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-comparative-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
