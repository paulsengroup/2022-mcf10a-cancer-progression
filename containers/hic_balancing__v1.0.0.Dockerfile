# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM curlimages/curl:7.87.0 AS downloader

ARG CONTAINER_VERSION
ARG JUICERTOOLS_VER=2.20.00
ARG HICTOOLS_VER=3.30.00

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi


ARG JUICERTOOLS_URL="https://github.com/aidenlab/Juicebox/releases/download/v${JUICERTOOLS_VER}/juicer_tools.${JUICERTOOLS_VER}.jar"
ARG HICTOOLS_URL="https://github.com/aidenlab/HiCTools/releases/download/v${HICTOOLS_VER}/hic_tools.${HICTOOLS_VER}.jar"

ARG JUICERTOOLS_SHA256='5de743337191a1b6b5bbc41daae4599766d4bf7073410ca2e1aa29376b306545'
ARG HICTOOLS_SHA256='2b09b0642a826ca5730fde74e022461a708caf62ed292bc5baaa841946721867'

RUN cd /tmp \
&& curl -LO "$JUICERTOOLS_URL" \
&& curl -LO "$HICTOOLS_URL" \
&& curl -L 'https://raw.githubusercontent.com/aidenlab/juicer/1c414ddebc827849f9c09fe2d2a2ea7c9a8c78df/LICENSE' -o juicer_tools.LICENSE \
&& curl -L 'https://raw.githubusercontent.com/aidenlab/HiCTools/6b2fab8e78685deae199c33bbb167dcab1dbfbb3/LICENSE' -o hic_tools.LICENSE \
&& echo "$JUICERTOOLS_SHA256  $(basename "$JUICERTOOLS_URL")" > checksum.sha256 \
&& echo "$HICTOOLS_SHA256  $(basename "$HICTOOLS_URL")" >> checksum.sha256 \
&& sha256sum -c checksum.sha256 \
&& chmod 644 *.jar *LICENSE


FROM ghcr.io/paulsengroup/ci-docker-images/ubuntu-22.04-cxx-gcc-12:20230204 AS builder

COPY containers/assets/hic2cool-ng-de0fad0.tar.xz /tmp/

RUN apt-get update \
&& apt-get install -y libtbb2-dev
RUN cd /tmp \
&& tar -xf hic2cool-ng-*.tar.xz \
&& cmake -DCMAKE_BUILD_TYPE=Release \
         -S /tmp/hic2cool-ng*/ \
         -B /tmp/build \
&& cmake --build /tmp/build -j "$(nproc)"

FROM ubuntu:22.04 AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION
ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y git \
                      libcurl4 \
                      libcurl4-openssl-dev \
                      libtbb2 \
                      openjdk-18-jre \
                      pigz \
                      python3 \
                      python3-pip \
                      zstd \
&& pip install 'bioframe>=0.3.3' \
                git+https://github.com/robomics/cooler.git@balance-cis-bugfix \
               'hic2cool>=0.8.3' \
               'hic-straw>=1.3.1' \
&& apt-get remove -y git \
                     libcurl4-openssl-dev \
                     python3-pip \
&& rm -rf /var/lib/apt/lists/*

COPY --from=downloader --chown=root:root /tmp/juicer_tools*.jar /usr/local/share/java/juicer_tools/
COPY --from=downloader --chown=root:root /tmp/hic_tools*.jar /usr/local/share/java/hic_tools/
COPY --from=builder --chown=root:root /tmp/build/src/hic2cool-ng /usr/local/bin/hic2cool-ng
COPY --from=builder --chown=root:root /tmp/hic2cool-ng*/utils/cool2hic.py /usr/local/bin/cool2hic-ng

COPY --from=downloader --chown=root:root /tmp/juicer_tools.LICENSE /usr/local/share/licenses/juicer_tools/LICENSE
COPY --from=downloader --chown=root:root /tmp/hic_tools.LICENSE /usr/local/share/licenses/hic_tools/LICENSE
COPY --from=builder --chown=root:root /tmp/hic2cool-ng*/LICENSE /usr/local/share/hic2cool-ng/LICENSE

RUN printf '%s\nexec /usr/bin/java -Xms512m -Xmx16g -jar %s "$@"\n' \
      '#!/bin/sh' \
      "$(printf '%s' /usr/local/share/java/juicer_tools/juicer_tools*.jar)" > /usr/local/bin/juicer_tools \
&& printf '%s\nexec /usr/bin/java -Xms512m -Xmx16g -jar %s "$@"\n' \
      '#!/bin/sh' \
      "$(printf '%s' /usr/local/share/java/hic_tools/hic_tools*.jar)" > /usr/local/bin/hic_tools \
&& chmod 755 /usr/local/bin/*_tools

CMD ["/bin/bash"]
WORKDIR /data

RUN cooler --help
RUN hic2cool --help
RUN hic2cool-ng --help
RUN cool2hic-ng --help
RUN juicer_tools --help
RUN hic_tools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hic_balancing}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
