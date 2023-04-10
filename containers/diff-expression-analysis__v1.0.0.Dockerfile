# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.4.2 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

USER root
RUN apt-get update \
&& apt-get install -y less \
&& rm -rf /var/lib/apt/lists/*
USER mambauser

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DESEQ2_VERSION="1.38.*"
ARG CLUSTERPROFILER_VERSION="4.6.*"
ARG ENRICHPLOT_VERSION="1.18.*"
ARG ORG_HS_EG_DB="3.16.*"

RUN micromamba install -y                                              \
               -c conda-forge                                          \
               -c bioconda                                             \
               "bioconductor-clusterprofiler=$CLUSTERPROFILER_VERSION" \
               "bioconductor-deseq2=$DESEQ2_VERSION"                   \
               "bioconductor-enrichplot=$ENRICHPLOT_VERSION"           \
               "bioconductor-org.hs.eg.db=$ORG_HS_EG_DB"               \
               bioframe                                                \
               gprofiler-official                                      \
               r-argparse                                              \
               r-fs                                                    \
               r-optparse                                              \
               r-pheatmap                                              \
               r-stringr                                               \
               procps-ng                                               \
&& micromamba clean --all -y

RUN touch /opt/conda/lib/R/etc/.Rprofile

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
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

RUN Rscript --no-save -e 'quit(status=!library("clusterProfiler", character.only=T, logical.return=T), save="no")'
RUN Rscript --no-save -e 'quit(status=!library("DESeq2", character.only=T, logical.return=T), save="no")'
RUN Rscript --no-save -e 'quit(status=!library("enrichplot", character.only=T, logical.return=T), save="no")'
RUN Rscript --no-save -e 'quit(status=!library("org.Hs.eg.db", character.only=T, logical.return=T), save="no")'

RUN python3 -c 'import bioframe; import gprofiler'

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-diff-expression-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
