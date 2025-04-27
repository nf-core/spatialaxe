FROM mambaorg/micromamba:1.5.10-noble
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/conda.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y
# RUN echo "channel_alias: https://conda-mirror.artifacts.altoslabs.com/" >> $HOME/.condarc


RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
