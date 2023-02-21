# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS downloader
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_NO_CACHE_DIR=0

ARG CONTAINER_VERSION
ARG DCHIC_VER=${CONTAINER_VERSION}

COPY "containers/assets/dchic-${DCHIC_VER}.tar.xz" /tmp/

RUN apt-get update \
&& apt-get install -y findutils sed tar xz-utils \
&& cd /tmp \
&& tar -xf dchic-*.tar.xz \
&& mv dchic-*/ dchic \
&& rm -r dchic/demo dchic/docs \
&& sed -i 's/name: dchic/name: base/' dchic/packages/dchic.yml \
&& find dchic -type f -exec chmod uga+r {} + \
&& find dchic -type d -exec chmod uga+rx {} + \
&& chmod 755 dchic/dchicf.r

FROM mambaorg/micromamba:1.3.1 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DCHIC_VERSION="${CONTAINER_VERSION}"

COPY --from=downloader --chown=nobody:nogroup /tmp/dchic /opt/dchic

RUN cd /opt/dchic \
&& micromamba install -y -f packages/dchic.yml \
&& micromamba clean --all -y \
&& R CMD INSTALL packages/functionsdchic_*.tar.gz

# Workaround for ImportError raised by from pysam.libchtslib import * from igv_reports/datauri.py
# ImportError: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
# I am pretty sure we're not actually accessing any symbol from libcrypto, so this workaround should be ok.
RUN cd /opt/conda/lib && ln -s libcrypto.so libcrypto.so.1.0.0

RUN touch /opt/conda/lib/R/etc/.Rprofile

USER root
RUN mkdir /opt/dchic/bin \
&& ln -s /opt/dchic/dchicf.r /opt/dchic/bin/dchicf.r \
&& chown -R nobody:nogroup /opt/dchic
USER mambauser


WORKDIR /data

ENV PATH="/opt/conda/bin:/opt/dchic/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/dchic/bin/dchicf.r"]
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

RUN dchicf.r --help
RUN fithic --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-dchic}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
