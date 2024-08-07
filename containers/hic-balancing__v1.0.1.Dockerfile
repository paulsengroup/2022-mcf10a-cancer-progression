# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM curlimages/curl:8.1.0 AS downloader

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


FROM ghcr.io/paulsengroup/ci-docker-images/ubuntu-22.04-cxx-clang-15:20230501 AS hic2cool-ng

RUN apt-get update \
&& apt-get install -y libtbb2-dev

COPY containers/assets/hic2cool-ng-6e27c65.tar.xz /tmp/


RUN cd /tmp \
&& tar -xf hic2cool-ng-*.tar.xz                        \
&& cd hic2cool-ng*/                                    \
&& sed -i 's|hdf5/.*|hdf5/1.14.3|' conanfile.txt        \
&& sed -i 's|highfive/.*|highfive/2.9.0|' conanfile.txt \
&& sed -i 's|zlib/.*|zlib/1.3.1|' conanfile.txt         \
&& conan install .                                     \
    --build=missing                                    \
    --build=cascade                                    \
    -pr:b="$CONAN_DEFAULT_PROFILE_PATH"                \
    -pr:h="$CONAN_DEFAULT_PROFILE_PATH"                \
    -s build_type=Release                              \
    -s compiler.cppstd=20                              \
    --output-folder=build

RUN cd /tmp/hic2cool-ng*/ \
&& cmake -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_PREFIX_PATH="$PWD/build" \
         -DCMAKE_INSTALL_PREFIX=/tmp/hic2cool-ng \
         -S . \
         -B build/ \
&& cmake --build build -j $(nproc) \
&& cmake --install build


FROM mambaorg/micromamba:1.5.8-jammy AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION
ARG PIP_NO_CACHE_DIR=0
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# https://github.com/open2c/cooler/pull/313
# https://github.com/open2c/cooler/pull/323
# https://github.com/open2c/cooler/pull/324
COPY containers/patches/cooler.patch /tmp/

RUN micromamba install \
      -c conda-forge \
      -c bioconda \
      'bioframe>=0.4.1' \
      'hic2cool>=0.8.3' \
      'hic-straw>=1.3.1' \
      'hictkpy==0.0.5' \
      'iced==0.5.10' \
      'numpy<2' \
&& micromamba clean --all -y

USER root
RUN apt-get update \
&& apt-get install -y git \
                      libcurl4 \
                      libcurl4-openssl-dev \
                      libtbb2 \
                      openjdk-18-jre \
                      pigz \
                      procps \
                      zstd \
&& git clone https://github.com/open2c/cooler.git /tmp/cooler \
&& cd /tmp/cooler \
&& git checkout v0.9.1 \
&& patch -p0 < /tmp/cooler.patch \
&& pip install --upgrade pip setuptools \
&& pip install /tmp/cooler \
&& apt-get remove -y git \
                     libcurl4-openssl-dev \
&& rm -rf /var/lib/apt/lists/* \
          /tmp/cooler* \

USER mambauser

COPY --from=downloader  --chown=root:root /tmp/juicer_tools*.jar              /usr/local/share/java/juicer_tools/
COPY --from=downloader  --chown=root:root /tmp/hic_tools*.jar                 /usr/local/share/java/hic_tools/
COPY --from=hic2cool-ng --chown=root:root /tmp/hic2cool-ng/bin/hic2cool-ng    /usr/local/bin/hic2cool-ng
COPY --from=hic2cool-ng --chown=root:root /tmp/hic2cool-ng*/utils/cool2hic.py /usr/local/bin/cool2hic-ng

COPY --from=downloader  --chown=root:root /tmp/juicer_tools.LICENSE           /usr/local/share/licenses/juicer_tools/LICENSE
COPY --from=downloader  --chown=root:root /tmp/hic_tools.LICENSE              /usr/local/share/licenses/hic_tools/LICENSE
COPY --from=hic2cool-ng --chown=root:root /tmp/hic2cool-ng*/LICENSE           /usr/local/share/hic2cool-ng/LICENSE

RUN printf '%s\nexec /usr/bin/java -Xms512m -Xmx16g -jar %s "$@"\n' \
      '#!/bin/sh' \
      "$(printf '%s' /usr/local/share/java/juicer_tools/juicer_tools*.jar)" > /usr/local/bin/juicer_tools \
&& printf '%s\nexec /usr/bin/java -Xms512m -Xmx16g -jar %s "$@"\n' \
      '#!/bin/sh' \
      "$(printf '%s' /usr/local/share/java/hic_tools/hic_tools*.jar)" > /usr/local/bin/hic_tools \
&& chmod 755 /usr/local/bin/*_tools

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

RUN cooler --help
RUN hic2cool --help
RUN hic2cool-ng --help
RUN cool2hic-ng --help
RUN juicer_tools --help
RUN hic_tools --help

RUN python3 -c 'import bioframe; import hictkpy; import hicstraw; import iced'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hic_balancing}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
