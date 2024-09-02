# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.5.9 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y \
               -c conda-forge \
               h5py \
               matplotlib \
               'numpy>=2' \
               'opencv=4.10.*' \
               'pandas>=2' \
               'pillow=10.3.*' \
               'scipy=1.14.*' \
               'scikit-image=0.24.*' \
               procps-ng \
&& micromamba clean --all -y

USER root
RUN apt-get update \
&& apt-get install -y libegl1 libgl1 libopengl0 \
&& rm -rf /var/lib/apt/lists/*
USER mambauser

ENV MKL_NUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1
ENV OMP_NUM_THREADS=1

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

RUN python3 -c 'import h5py, matplotlib, numpy, cv2, pandas, PIL, scipy, skimage'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-image-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
