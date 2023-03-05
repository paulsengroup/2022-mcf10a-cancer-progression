# Copyright (C) 2023 Saleh Oshaghi <mohao@uio.no>
#
# SPDX-License-Identifier: MIT
FROM ubuntu:22.04 AS builder

#ARG Chrom3D_URL="https://github.com/Chrom3D/Chrom3D/archive/refs/tags/v1.0.2.tar.gz"
#ARG Chrom3D_CHECKSUM="f19e046216f17c43369ccddf45cdc114b364cf3ad7173501a56ab87a7d310be9  v1.0.2.tar.gz"


RUN apt-get update \
&& apt-get install -y cmake \
                      g++ \
                      git \
                      libboost-dev \
                      libboost-filesystem-dev \
                      libboost-random-dev

RUN git config --global user.email bot@example.com \
&& git config --global user.name bot

RUN git clone 'https://github.com/Chrom3D/Chrom3D.git' /tmp/Chrom3D \
&& cd /tmp/Chrom3D \
&& git remote add robomics 'https://github.com/robomics/Chrom3D.git' \
&& git fetch robomics \
&& git switch -c patched v1.0.2 \
&& git merge robomics/add-cmake-support -m 'patch #1' \
&& git merge robomics/fix-compiler-warnings -m 'patch #2'

RUN cmake -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX=/tmp/chrom3d-staging \
         -S /tmp/Chrom3D/ \
         -B /tmp/build/ \
&& cmake --build /tmp/build/ -j $(nproc) \
&& cmake --install /tmp/build/

#&& echo "$Chrom3D_CHECKSUM" >> checksums.sha256 \
#&& shasum -a256 -c checksums.sha256 \
#&& tar -xf "v1.0.2.tar.gz"


FROM ubuntu:22.04 AS base

RUN ln -snf /usr/share/zoneinfo/CET /etc/localtime \
&& echo CET | tee /etc/timezone > /dev/null

ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y bedtools \
                      gawk \
                      python3 \
                      python3-pip \
&& pip install 'pandas==1.5.*' \
&& apt-get remove -y python3-pip \
&& apt-get autoremove -y \
&& rm -rf /var/lib/apt/lists/*


ARG CONTAINER_VERSION
ARG CONTAINER_TITLE


COPY --from=builder "/tmp/chrom3d-staging/" "/usr/local/"

CMD ["/usr/local/bin/Chrom3D"]
WORKDIR /data

RUN bedtools --version
RUN Chrom3D --version
RUN python3 -c 'import pandas; print(pandas.__version__)'


LABEL org.opencontainers.image.authors='Saleh Oshaghi <mohao@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-Chrom3D}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
