# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS builder

ARG CONTAINER_VERSION
ARG SAMTOOLS_VER=${CONTAINER_VERSION}
ARG SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2"
ARG SAMTOOLS_CHECKSUM="2fa0a25f78594cf23d07c9d32d5060a14f1c5ee14d7b0af7a8a71abc9fdf1d07  samtools-$SAMTOOLS_VER.tar.bz2"

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
&& apt-get install -y autoconf \
                      bzip2 \
                      curl \
                      gawk \
                      g++ \
                      make \
                      perl \
                      libbz2-dev \
                      liblzma-dev \
                      zlib1g-dev

RUN cd /tmp \
&& curl -LO "$SAMTOOLS_URL" \
&& echo "$SAMTOOLS_CHECKSUM" >> checksums.sha256 \
&& shasum -a256 -c checksums.sha256 \
&& tar -xf "samtools-$SAMTOOLS_VER.tar.bz2"

RUN cd "/tmp/samtools-$SAMTOOLS_VER" \
&& LDFLAGS='-static' ./configure --enable-configure-htslib \
                                 --without-curses          \
                                 --prefix=/tmp/staging     \
&& make -j "$(nproc)" \
&& make -j "$(nproc)" test \
&& make install \
&& install -Dm0644 LICENSE /tmp/staging/share/doc/samtools/copyright

FROM ubuntu:22.04 AS base
ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

COPY --from=builder "/tmp/staging/bin" "/usr/local/bin"
COPY --from=builder "/tmp/staging/share" "/usr/local/share"

ENTRYPOINT ["/usr/local/bin/samtools"]
WORKDIR /data

RUN samtools --version
RUN samtools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-samtools}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
