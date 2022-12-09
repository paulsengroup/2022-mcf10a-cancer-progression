# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS builder

ARG CONTAINER_VERSION
ARG BEDTOOLS_VER=${CONTAINER_VERSION}
ARG BEDTOOLS_URL="https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VER/bedtools-$BEDTOOLS_VER.tar.gz"

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
&&  apt-get install -y curl \
                       diffutils \
                       g++ \
                       make \
                       python2 \
                       libbz2-dev \
                       liblzma-dev \
                       zlib1g-dev

RUN cd /tmp && curl -L "$BEDTOOLS_URL" | tar -xzf -

RUN mkdir -p /tmp/bin \
&&  ln -s /usr/bin/python2 /tmp/bin/python

RUN cd "/tmp/bedtools2" \
&&  PATH="/tmp/bin:$PATH" make -j $(nproc)  \
&&  PATH="/tmp/bin:$PATH" make test 2>&1 | grep -q 'Tools failing:  negativecontrol'

FROM ubuntu:22.04 AS base
ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

COPY --from=builder "/tmp/bedtools2/bin" "/usr/local/bin"
COPY --from=builder "/tmp/bedtools2/LICENSE" "/usr/local/share/licenses/bedtools/LICENSE"

RUN apt-get update \
&&  apt-get install -y bzip2 \
                       xz-utils \
                       zlib1g \
&&  rm -rf /var/lib/apt/lists/*

RUN chown -R root:root /usr/local/share/licenses/

CMD ["/usr/local/bin/bedtools"]
WORKDIR /data

RUN bedtools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-bedtools}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
