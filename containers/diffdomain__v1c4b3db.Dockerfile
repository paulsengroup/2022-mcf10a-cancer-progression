# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
FROM ubuntu:22.04 AS downloader

ARG CONTAINER_VERSION
ARG DIFFDOMAIN_VER=${CONTAINER_VERSION}

COPY "containers/assets/robomics-diffdomain-${DIFFDOMAIN_VER}.tar.xz" /tmp/

RUN apt-get update \
&& apt-get install -y tar xz-utils \
&& cd /tmp \
&& tar -xf robomics-diffdomain-*.tar.xz \
&& mv robomics-diffdomain-*/ diffdomain \
&& chmod 755 diffdomain/diffdomain-py3/diffdomains.py \
             diffdomain/diffdomain-py3/classification.py \
&& rm diffdomain/diffdomain-py3/__init__.py \
&& find diffdomain -type f -exec chmod uga+r {} + \
&& find diffdomain -type d -exec chmod uga+rx {} +

FROM mambaorg/micromamba:2.0.3 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG DIFFDOMAIN_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        'python=3.10' \
        'hic-straw=1.3.1' \
        'cooler>=0.9' \
        docopt \
        h5py \
        matplotlib \
        numpy \
        'pandas<2' \
        procps-ng \
        seaborn \
        statsmodels \
        tqdm \
        tracywidom \
&& micromamba clean --all -y

COPY --from=downloader --chown=nobody:nogroup \
    /tmp/diffdomain/diffdomain-py3/*.py \
    /opt/diffdomain_py3/lib/

COPY --from=downloader --chown=nobody:nogroup \
    /tmp/diffdomain/LICENSE \
    /opt/diffdomain_py3/share/

USER root
RUN mkdir /opt/diffdomain_py3/bin \
&& ln -s /opt/diffdomain_py3/lib/diffdomains.py /opt/diffdomain_py3/bin/diffdomains.py \
&& ln -s /opt/diffdomain_py3/lib/classification.py /opt/diffdomain_py3/bin/classification.py \
&& chown -R nobody:nogroup /opt/diffdomain_py3
USER mambauser

ENV PATH="/opt/conda/bin:/opt/diffdomain_py3/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/diffdomain_py3/bin/diffdomains.py"]
WORKDIR /data


RUN diffdomains.py --help
RUN classification.py --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-diffdomains}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
