# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM curlimages/curl:7.87.0 AS downloader

ARG CONTAINER_VERSION
ARG JUICER_VER=1.22.01

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG JUICER_URL="https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_${JUICER_VER}.jar"
ARG JUICER_SHA256="5bd863e1fbc4573de09469e0adc5ab586e2b75b14dd718465e14dc299d7243a0"

RUN cd /tmp \
&& curl -LO "$JUICER_URL" \
&& curl -L 'https://raw.githubusercontent.com/aidenlab/juicer/1c414ddebc827849f9c09fe2d2a2ea7c9a8c78df/LICENSE' -o LICENSE \
&& echo "$JUICER_SHA256  $(basename "$JUICER_URL")" > checksum.sha256 \
&& sha256sum -c checksum.sha256 \
&& chmod 644 juicer_tools*.jar LICENSE


FROM ubuntu:22.04 AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION
ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&&  apt-get install -y cython3 \
                       openjdk-18-jre \
                       pigz \
                       python3 \
                       python3-pip \
                       zstd \
&& pip install "matplotlib<3.6" \
               "numpy<1.24" \
                cooler \
&& apt-get remove -y cython3 python3-pip \
&& rm -rf /var/lib/apt/lists/*

COPY --from=downloader --chown=root:root /tmp/juicer_tools*.jar /usr/local/share/java/juicer_tools/
COPY --from=downloader --chown=root:root /tmp/LICENSE /usr/local/share/licenses/juicer_tools/LICENSE

RUN printf '%s\nexec /usr/bin/java -Xms512m -Xmx2048m -jar %s "$@"\n' \
      '#!/bin/sh' \
      "$(printf '%s' /usr/local/share/java/juicer_tools/juicer_tools*.jar)" > /usr/local/bin/juicer_tools \
&& chmod 755 /usr/local/bin/juicer_tools

CMD ["/usr/local/bin/juicer_tools"]
WORKDIR /data

RUN cooler --help
RUN juicer_tools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-cooler_to_hic}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
