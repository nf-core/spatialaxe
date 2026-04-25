
# Dockerfile to create container for SCS
# SCS: cell segmentation on spatial transcriptomics + morphology

FROM condaforge/mambaforge:latest

LABEL authors="Katarzyna Wreczycka" \
    description="Image for SCS_SEGMENT module"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    git \
    curl \
    libhdf5-dev \
    pkg-config \
    libgl1 \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

ARG SCS_COMMIT=eccc45049c7d5ca4fb9717abe913af75e3325173
RUN git clone https://github.com/katwre/SCS.git /app/SCS && \
    cd /app/SCS && \
    git checkout ${SCS_COMMIT}

RUN pip config set global.timeout 300 && pip config set global.retries 20

RUN mamba env create -f /app/SCS/environment.yml && \
    mamba clean --all -f -y

# Ensure all pip installs are in the SCS conda environment
SHELL ["conda", "run", "-n", "SCS", "/bin/bash", "-c"]
RUN pip install --no-cache-dir spateo-release==1.1.0 scikit-misc==0.3.1
RUN pip install --no-cache-dir pandas numpy tifffile pyarrow
RUN pip install --no-cache-dir alphashape

## fastpd install removed: repo does not exist and is not required

ENV CONDA_DEFAULT_ENV=SCS \
    PATH=/opt/conda/envs/SCS/bin:$PATH \
    PYTHONPATH=/app/SCS

# drop directly into the SCS env shell 
CMD ["bash", "--noprofile", "--norc", "-lc", "source /opt/conda/etc/profile.d/conda.sh && conda activate SCS && export PS1=\"(SCS) \\u@\\h:\\w\\$ \" && exec bash --noprofile --norc -i"]


