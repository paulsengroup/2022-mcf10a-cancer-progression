# Copyright (C) 2023 Saleh Oshaghi <mohao@uio.no>
#
# SPDX-License-Identifier: MIT
FROM ubuntu:22.04 AS builder

#ARG Chrom3D_URL="https://github.com/Chrom3D/Chrom3D/archive/refs/tags/v1.0.2.tar.gz"
#ARG Chrom3D_CHECKSUM="f19e046216f17c43369ccddf45cdc114b364cf3ad7173501a56ab87a7d310be9  v1.0.2.tar.gz"


RUN apt-get update \
&& apt-get install -y autoconf \
                      bzip2 \
                      curl \
                      gawk \
                      g++ \
                      make \
                      python3 \
                      bedtools \
                      libbz2-dev \
                      liblzma-dev \
                      zlib1g-dev \
                      git \
                      libboost-dev


RUN cd /tmp \
&& git clone https://github.com/robomics/Chrom3D.git -b fix-compiler-warnings
#&& echo "$Chrom3D_CHECKSUM" >> checksums.sha256 \
#&& shasum -a256 -c checksums.sha256 \
#&& tar -xf "v1.0.2.tar.gz"

ARG SRCPATH=tmp/Chrom3D/src
ARG TCLAPINCLUDE=tmp/Chrom3D/src/tclap-1.2.1/include
ARG BOOSTINCLUDE=$(whereis boost)
ARG CFLAGS=-O3

RUN g++ -O3 -std=c++11 -I $TCLAPINCLUDE -I BOOSTINCLUDE \
     $SRCPATH/Util.cpp $SRCPATH/Bead.cpp $SRCPATH/Chromosome.cpp $SRCPATH/Randomizer.cpp \
     $SRCPATH/Constraint.cpp $SRCPATH/Model.cpp $SRCPATH/MCMC.cpp $SRCPATH/Chrom3D.cpp -o Chrom3D



FROM ubuntu:22.04 AS base

RUN ln -snf /usr/share/zoneinfo/CET /etc/localtime \
&& echo CET | tee /etc/timezone > /dev/null

ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y  \
                      gawk \
                      python3 \
                      bedtools \
                      libboost-dev \
                      python3-dev             \
                      python3-pip             \
&& pip install pandas                            \
&&   rm -rf /var/lib/apt/lists/*



ARG CONTAINER_VERSION
ARG CONTAINER_TITLE


COPY --from=builder "Chrom3D" "/usr/local/bin"
CMD ["/usr/local/bin/Chrom3D"]
WORKDIR /data

RUN Chrom3D --help
RUN python3 --version


LABEL org.opencontainers.image.authors='Saleh Oshaghi <mohao@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-Chrom3D}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"