# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.3.1 AS base

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
        git \
        h5py \
        matplotlib \
        numpy \
        pandas \
        seaborn \
        statsmodels \
        tqdm \
        tracywidom \
&& git clone https://github.com/robomics/diffDomain.git /tmp/diffDomain \
&& cd /tmp/diffDomain \
&& git checkout "$DIFFDOMAIN_VER" \
&& micromamba remove -y git \
&& micromamba clean --all -y \
&& install -Dm0755 diffdomain-py3/diffdomains.py /opt/conda/lib/python3.10/site-packages/diffdomain_py3/diffdomains.py \
&& install -Dm0755 diffdomain-py3/classification.py /opt/conda/lib/python3.10/site-packages/diffdomain_py3/classification.py \
&& install -Dm0644 diffdomain-py3/utils.py /opt/conda/lib/python3.10/site-packages/diffdomain_py3/utils.py \
&& ln -s /opt/conda/lib/python3.10/site-packages/diffdomain_py3/diffdomains.py /opt/conda/bin/diffdomains.py \
&& ln -s /opt/conda/lib/python3.10/site-packages/diffdomain_py3/classification.py /opt/conda/bin/classification.py \
&& cd / \
&& rm -r /tmp/diffDomain

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/opt/conda/bin/diffdomains"]
WORKDIR /data


RUN diffdomains.py --help
RUN classification.py --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-diffdomains}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

