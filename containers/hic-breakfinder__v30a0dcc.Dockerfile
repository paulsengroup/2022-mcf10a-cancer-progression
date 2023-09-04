# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.5.0 AS builder

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG HIC_BREAKFINDER_GIT='https://github.com/dixonlab/hic_breakfinder.git'
ARG HIC_BREAKFINDER_TAG="${CONTAINER_VERSION}"

RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    bamtools \
    eigen \
    git \
    gxx \
    make

RUN cd /tmp \
&& git clone "$HIC_BREAKFINDER_GIT" \
&& cd hic_breakfinder \
&& git checkout "$HIC_BREAKFINDER_TAG"

RUN cd hic_breakfinder \
&& CXXFLAGS='-isystem /opt/conda/include/bamtools -isystem /opt/conda/include/eigen3 -fpermissive' \
  ./configure \
&& make -j $(nproc)


FROM mambaorg/micromamba:1.5.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG HINT_VERSION="${CONTAINER_VERSION}"


RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    bamtools \
    procps-ng \
&& micromamba clean --all -y

COPY --from=builder --chown=nobody:nogroup /tmp/hic_breakfinder/src/hic_breakfinder /usr/local/bin/
COPY --from=builder --chown=nobody:nogroup /tmp/hic_breakfinder/LICENSE /usr/local/share/licenses/hic_breakfinder/

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/usr/local/bin/hic_breakfinder"]
WORKDIR /data

RUN whereis hic_breakfinder

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hic_breakfinder}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
