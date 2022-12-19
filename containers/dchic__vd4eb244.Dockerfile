# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:1.1.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DC_HIC_VERSION="${CONTAINER_VERSION}"

ARG DC_HIC_GIT="https://github.com/ay-lab/dcHiC.git"
ARG DC_HIC_GIT_SHA="${CONTAINER_VERSION}"

RUN micromamba install -y             \
               -c conda-forge         \
               git                    \
&& git clone "$DC_HIC_GIT" /tmp/dchic \
&& cd /tmp/dchic                      \
&& git checkout "$DC_HIC_GIT_SHA"     \
&& sed -i 's/name: dchic/name: base/' packages/dchic.yml                \
&& micromamba install -y -f packages/dchic.yml                          \
&& micromamba remove -y git                                             \
&& micromamba clean --all -y                                            \
&& R CMD INSTALL packages/functionsdchic_*.tar.gz                       \
&& install -Dm0755 dchicf.r /opt/conda/bin/dchicf.r                     \
&& install -Dm0755 utility/*.{r,py,sh} /opt/conda/bin/                  \
&& install -Dm0755 utility/Chromosome_ArmWise_PCA/*.pl /opt/conda/bin/  \
&& install -Dm0644 LICENSE /opt/conda/share/licenses/dchic/LICENSE      \
&& cd / && rm -rf /tmp/dchic

# Workaround for ImportError raised by from pysam.libchtslib import * from igv_reports/datauri.py
# ImportError: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
# I am pretty sure we're not actually accessing any symbol from libcrypto, so this workaround should be ok.
RUN cd /opt/conda/lib && ln -s libcrypto.so libcrypto.so.1.0.0

RUN touch /opt/conda/lib/R/etc/.Rprofile

USER root
RUN chown -R nobody:nogroup /opt/conda/bin  /opt/conda/share/licenses/dchic/
USER mambauser


WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/usr/local/bin/dchicf.r"]
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
