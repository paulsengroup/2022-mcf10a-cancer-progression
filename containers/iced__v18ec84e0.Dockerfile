# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.2.0 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG ICED_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        'python>3' \
        bioframe \
        c-compiler \
        cooler \
        cython \
        git \
        make \
        'matplotlib<3.6' \
        'numpy>1.16,<1.24' \
        pandas \
        scikit-learn \
        'scipy>0.19' \
&& git clone https://github.com/robomics/iced.git /tmp/iced \
&& cd /tmp/iced \
&& git checkout "$ICED_VER" \
&& make cython -j "$(nproc)" \
&& pip install . \
&& cd / \
&& rm -r /tmp/iced \
&& micromamba remove -y \
        c-compiler \
        cython \
        git \
        make \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data


RUN cooler --help
RUN python3 -c "import cooler; import iced"

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-iced}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

