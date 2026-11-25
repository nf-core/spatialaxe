# Build from the repository root:
# docker build -t split:0.3.0 modules/local/split
FROM mambaorg/micromamba:2.3.3@sha256:800e7ade3ffe29c9a9ac2026163131495f8197c3852e572c5835beb4e8a33cd6

USER root
RUN apt-get update && apt-get install --yes --no-install-recommends procps gzip && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY build-environment.yml /tmp/build-environment.yml
RUN micromamba create --yes --strict-channel-priority --file /tmp/build-environment.yml && micromamba clean --all --yes

ENV ENV_NAME=split
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript -e 'remotes::install_github("bdsc-tds/SPLIT@e880e39c03a7036e11c9867701c76fdb26acc14c", dependencies=FALSE, upgrade="never", build_vignettes=FALSE); library(SPLIT); stopifnot(packageVersion("SPLIT") == "0.3.0")'

COPY --chmod=755 purify/resources/purify.R /usr/local/bin/split-purify
# Nextflow overrides the image entrypoint, so activation cannot depend on it.
ENV PATH="/opt/conda/envs/split/bin:${PATH}"
