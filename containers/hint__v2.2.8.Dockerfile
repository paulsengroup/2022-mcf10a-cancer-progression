# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS download

ARG BICSEQ2_SEG_URL='http://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2/BICseq2-seg_v0.7.3.tar.gz'

RUN apt-get update \
&& apt-get install -y curl \
&& cd /tmp \
&& curl -L "$BICSEQ2_SEG_URL" | tar -xzf -


FROM mambaorg/micromamba:2.0.8 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG HINT_VERSION="${CONTAINER_VERSION}"


RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    "hint=$HINT_VERSION" \
    bioframe \
    biopython \
    'pandas<1' \
    pysam \
    procps-ng \
    unzip \
&& micromamba clean --all -y

COPY --from=download --chown=nobody:nogroup /tmp/BICseq2-seg_v0.7.3 /opt/BICseq2-seg/
COPY --chown=nobody:nogroup containers/patches/hint.patch /tmp/

USER root
RUN apt-get update \
&& apt-get install -y patch \
&& patch /opt/conda/bin/hint -i /tmp/hint.patch \
&& rm /tmp/hint.patch \
&& apt-get remove -y patch \
&& apt-get autoremove -y \
&& rm -rf /var/lib/apt/lists/*
USER mambauser

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/hictrans/bin/hint"]
WORKDIR /data


# We have to explicitly set these R_* env variables in order for the
# container to work correctly when running using Apptainer
ENV R_HOME=/opt/conda/lib/R
ENV R_LIBS=/opt/conda/lib/R/lib
ENV R_ENVIRON=/opt/conda/lib/R/etc/Renviron
ENV R_HISTFILE=/tmp/.Rhistory

ENV R_HOME_USER='$R_HOME'
ENV R_LIBS_USER='$R_LIBS'
ENV R_ENVIRON_USER='$R_ENVIRON'
ENV R_PROFILE_USER=/opt/conda/lib/R/etc/.Rprofile

RUN hint --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hint}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
