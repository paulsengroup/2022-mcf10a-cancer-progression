# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.0.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DC_HIC_VERSION="${CONTAINER_VERSION}"

ARG DC_HIC_GIT="https://github.com/robomics/dcHiC.git"
ARG DC_HIC_GIT_SHA="${CONTAINER_VERSION}"

RUN micromamba install -y             \
               -c conda-forge         \
               git                    \
&& git clone "$DC_HIC_GIT" /tmp/dchic \
&& cd /tmp/dchic                      \
&& git checkout "$DC_HIC_GIT_SHA"     \
&& micromamba remove -y git           \
&& sed -i 's/name: dchic/name: base/' packages/dchic.yml                \
&& micromamba install -y -f packages/dchic.yml                          \
&& micromamba clean --all -y                                            \
&& R CMD INSTALL packages/functionsdchic_*.tar.gz                       \
&& install -Dm0755 dchicf.r /opt/conda/bin/dchicf.r                     \
&& install -Dm0755 utility/*.{r,py,sh} /opt/conda/bin/                  \
&& install -Dm0755 utility/Chromosome_ArmWise_PCA/*.pl /opt/conda/bin/  \
&& cd / && rm -rf /tmp/dchic

USER root
RUN chown nobody:nogroup /opt/conda/bin/*
USER mambauser

RUN touch /opt/conda/lib/R/etc/.Rprofile

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/usr/local/bin/dchicf.r"]

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
