# GPU image for segger segmentation
# Base: PyTorch with CUDA 12.4 (consistent CUDA support)
FROM pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    procps \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Clone and install segger from fork with CUDA 12 support
# Using fork with fixed API (is_token_based, num_node_features) and data module exports
WORKDIR /workspace
RUN git clone https://github.com/an-altosian/segger_dev.git /workspace/segger_dev && \
    cd /workspace/segger_dev && \
    pip install --no-cache-dir -e ".[cuda12]"

# Install missing dependencies for data.io module and training
RUN pip install --no-cache-dir dask-geopandas lightning pytorch_lightning

# Fix squidpy/anndata compatibility - upgrade squidpy to version compatible with anndata 0.10+
RUN pip install --no-cache-dir "squidpy>=1.4.0"

# Make NVRTC discoverable for CuPy (needs libnvrtc.so.12 for JIT kernel compilation)
# The nvidia-cuda-nvrtc-cu12 pip package installs it under site-packages.
# Use both ldconfig AND symlinks for maximum reliability (ldconfig may not survive
# Wave container augmentation on Seqera Platform).
RUN NVRTC_LIB=$(python -c "import nvidia.cuda_nvrtc, pathlib; print(pathlib.Path(nvidia.cuda_nvrtc.__file__).parent / 'lib')") && \
    echo "$NVRTC_LIB" > /etc/ld.so.conf.d/nvidia-nvrtc.conf && \
    ln -sf "$NVRTC_LIB"/libnvrtc* /usr/lib/ && \
    ldconfig && \
    echo "NVRTC registered at: $NVRTC_LIB" && \
    ldconfig -p | grep nvrtc

# Install faiss-cpu for vector search (faiss-gpu not available via pip)
# Note: Removed cuvs-cu12 and nvidia library force-reinstalls as they break PyTorch CUDA detection
RUN pip install --no-cache-dir faiss-cpu

# Verify PyTorch was built with CUDA support (doesn't require GPU at build time)
RUN python -c "import torch; print(f'PyTorch version: {torch.__version__}'); print(f'CUDA compiled: {torch.version.cuda}')"

# Set environment variables for segger
ENV PYTHONPATH=/workspace/segger_dev/src:$PYTHONPATH
ENV PATH="$PATH:/workspace/segger_dev/src/segger/cli/"
ENV CUPY_CACHE_DIR="/tmp/cupy_cache"

# Set default shell
CMD ["/bin/bash"]
