# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.5.3 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y \
               -c conda-forge \
               -c bioconda \
               'bioframe>=0.4.1' \
               'hdbscan==0.8.*' \
               matplotlib \
               'networkx==3.*' \
               numpy \
               pandas \
               pigz \
               procps-ng \
               pyBigWig \
               scipy \
               seaborn \
&& micromamba clean --all -y

ENV MKL_NUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1
ENV OMP_NUM_THREADS=1

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

RUN python3 -c 'import bioframe, hdbscan, matplotlib, networkx, numpy, pandas, pyBigWig, scipy, seaborn'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-comparative-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
