# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS downloader

ARG GO_FIGURE_VERSION='1.0.1'
ARG GO_FIGURE_URL="https://gitlab.com/evogenlab/GO-Figure/-/archive/v${GO_FIGURE_VERSION}/GO-Figure-v${GO_FIGURE_VERSION}.tar.gz"
ARG GO_FIGURE_SHA256="d0169b62189800e7fae88ee475c8c19ec2034a5298e121f584ed53d54f54b942  GO-Figure-v${GO_FIGURE_VERSION}.tar.gz"

ARG GO_RELEASE='2022-11-03'
ARG GO_OBO_URL="http://release.geneontology.org/${GO_RELEASE}/ontology/go.obo"
ARG GOA_GAF_URL="http://release.geneontology.org/${GO_RELEASE}/annotations/goa_human.gaf.gz"

ARG GO_OBO_SHA256='e351a37a24168f106fca574accad4d27f20622e2b718db698a7e38515f0dba09  go.obo'
ARG GOA_GAF_SHA256='01fd1d5e1603f5714c62feefa2f853cf50ed2be63da60c803ac00b262b0490d8  goa_human.gaf.gz'

RUN printf '%s\n' "$GO_FIGURE_SHA256" "$GO_OBO_SHA256" "$GOA_GAF_SHA256" > /tmp/checksum.sha256

RUN cd /tmp/ \
&& apt-get update \
&& apt-get install -y curl diffutils patch \
&& curl -LO "$GO_FIGURE_URL" \
&& curl -LO "$GO_OBO_URL" \
&& curl -LO "$GOA_GAF_URL" \
&& sha256sum -c checksum.sha256 \
&& tar -xf *.tar.gz \
&& gzip -d goa_human.gaf.gz \
&& rm -f *.tar.gz

COPY containers/patches/gofigure.py.patch "/tmp/GO-Figure-v${GO_FIGURE_VERSION}/"

RUN cd "/tmp/GO-Figure-v${GO_FIGURE_VERSION}/" \
&& patch -b -p 1 -i gofigure.py.patch

RUN install -Dm0755 "/tmp/GO-Figure-v${GO_FIGURE_VERSION}/gofigure.py" /tmp/gofigure/gofigure.py \
&& install -Dm0755 "/tmp/GO-Figure-v${GO_FIGURE_VERSION}/scripts/"*.py /tmp/gofigure/ \
&& install -Dm0644 "/tmp/GO-Figure-v${GO_FIGURE_VERSION}/LICENSE" /tmp/gofigure/LICENSE

FROM mambaorg/micromamba:1.1.0 AS base

ARG CONTAINER_VERSION

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DESEQ2_VERSION="1.38.*"

RUN micromamba install -y                            \
               -c conda-forge                        \
               -c bioconda                           \
               adjusttext                            \
               bioframe                              \
               gprofiler-official                    \
               numpy                                 \
               scikit-learn                          \
               scipy                                 \
               seaborn                               \
&& micromamba clean --all -y

COPY --from=downloader --chown=nobody:nogroup /tmp/gofigure/*.py /opt/conda/bin/
COPY --from=downloader --chown=nobody:nogroup /tmp/gofigure/LICENSE /opt/conda/share/licenses/GO-Figure/

WORKDIR /data

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]

# Generate GO-Figure data
FROM base AS dataprep

COPY --from=downloader --chown=nobody:nogroup /tmp/go.obo /opt/conda/share/gofigure/data/
COPY --from=downloader --chown=nobody:nogroup /tmp/goa_human.gaf /opt/conda/share/gofigure/data/

USER root
RUN cd /opt/conda/share/gofigure/data/ \
&& relations.py go.obo > relations_full.tsv \
&& ics.py relations_full.tsv goa_human.gaf go.obo > ic.tsv \
&& chown -R nobody:nogroup /opt/conda/share/gofigure/data/
USER mambauser

FROM base as final

RUN python3 -c 'import bioframe; import gprofiler'
RUN gofigure.py --help

COPY --from=dataprep --chown=nobody:nogroup /opt/conda/share/gofigure/data/* /opt/conda/share/gofigure/data/
RUN ls -lah /opt/conda/share/gofigure/data/{ic,relations_full}.tsv

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.documentation='https://github.com/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2022-mcf10a-cancer-progression'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-functional-analysis}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
