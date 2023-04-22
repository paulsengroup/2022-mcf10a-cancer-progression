# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.4.2 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

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

RUN mkdir -p /opt/conda/lib/R/etc/ \
&& touch /opt/conda/lib/R/etc/.Rprofile

USER root
RUN apt-get update \
&& apt-get install -y less \
&& rm -rf /var/lib/apt/lists/*
USER mambauser

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DESEQ2_VERSION="1.38.*"
ARG APEGLM_VERSION="1.20.*"
ARG ASHR_VERSION="2.2_54"

RUN micromamba install -y                            \
               -c conda-forge                        \
               -c bioconda                           \
               "bioconductor-apeglm=$APEGLM_VERSION" \
               "bioconductor-deseq2=$DESEQ2_VERSION" \
               bioframe                              \
               numpy                                 \
               pandas                                \
               procps-ng                             \
               "r-ashr=$ASHR_VERSION"                \
               rpy2                                  \
&& micromamba clean --all -y

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
WORKDIR /data

# See https://www.rdocumentation.org/packages/bigparallelr/versions/0.3.1/topics/assert_cores
RUN echo 'try(bigparallelr::set_blas_ncores(1), silent = TRUE)' > /opt/conda/lib/R/etc/.Rprofile

RUN Rscript --no-save -e 'quit(status=!library("DESeq2", character.only=T, logical.return=T), save="no")'
RUN python3 -c 'import bioframe, rpy2'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-diff-expression-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
