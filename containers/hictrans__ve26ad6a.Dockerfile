# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS download

ARG CONTAINER_VERSION
ARG HICTRANS_VERSION="${CONTAINER_VERSION}"

RUN apt-get update \
&& apt-get install -y git \
&& cd /tmp \
&& git clone https://github.com/ay-lab/HiCtrans.git \
&& cd HiCtrans \
&& git checkout "$HICTRANS_VERSION"


FROM mambaorg/micromamba:2.0.3 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG HICTRANS_VERSION="${CONTAINER_VERSION}"


RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    cooler \
    pigz \
    procps-ng \
    r-catools \
    r-changepoint \
    r-data.table \
    r-depmixs4 \
    r-deoptimr \
    r-hashmap \
    r-optparse \
    r-r.utils \
    r-rcpp \
&& micromamba clean --all -y

COPY --from=download --chown=nobody:nogroup /tmp/HiCtrans/hictrans.v3.R /opt/hictrans/

USER root
RUN mkdir /opt/hictrans/bin \
&& printf '%s\n%s "$@"' \
    '#!/usr/bin/env sh' \
    'Rscript /opt/hictrans/hictrans.v3.R' > /opt/hictrans/bin/hictrans \
&& chown -R nobody:nogroup /opt/hictrans \
&& chmod 755 /opt/hictrans/bin/hictrans
USER mambauser

WORKDIR /data

ENV PATH="/opt/conda/bin:/opt/hictrans/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/hictrans/bin/hictrans"]
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

RUN hictrans --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hictrans}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
