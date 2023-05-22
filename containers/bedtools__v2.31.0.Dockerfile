# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS builder

ARG CONTAINER_VERSION
ARG BEDTOOLS_VER=${CONTAINER_VERSION}
ARG BEDTOOLS_URL="https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VER/bedtools-$BEDTOOLS_VER.tar.gz"
ARG BEDTOOLS_SHA256='8462980301da669b2976ef821a5cfd902c6653a5f7ed348a243a978ecf4593e1'
ARG BEDTOOLS_CHECKSUM="$BEDTOOLS_SHA256  bedtools-$BEDTOOLS_VER.tar.gz"

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
&&  apt-get install -y curl \
                       diffutils \
                       g++ \
                       make \
                       python3 \
                       libbz2-dev \
                       liblzma-dev \
                       zlib1g-dev

RUN cd /tmp \
&& curl -LO "$BEDTOOLS_URL" \
&& echo "$BEDTOOLS_CHECKSUM" > checksum.sha256 \
&& sha256sum -c checksum.sha256 \
&& tar -xzf *.tar.gz

RUN mkdir -p /tmp/bin \
&& ln -s /usr/bin/python3 /tmp/bin/python

RUN cd "/tmp/bedtools2" \
&& PATH="/tmp/bin:$PATH" make -j $(nproc)

RUN cd "/tmp/bedtools2" \
&& PATH="/tmp/bin:$PATH" make test 2>&1 | grep -q 'Tools failing:  negativecontrol'


FROM ubuntu:22.04 AS base
ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

COPY --from=builder "/tmp/bedtools2/bin" "/usr/local/bin"
COPY --from=builder "/tmp/bedtools2/LICENSE" "/usr/local/share/licenses/bedtools/LICENSE"

RUN apt-get update \
&&  apt-get install -y bzip2 \
                       procps \
                       xz-utils \
                       zlib1g \
&&  rm -rf /var/lib/apt/lists/*

RUN chown -R root:root /usr/local/share/licenses/

CMD ["/usr/local/bin/bedtools"]
WORKDIR /data

RUN bedtools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-bedtools}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
