# Dockerfile to create container with segtraq
# SegTraQ: A Python toolkit for quantitative and visual quality control
# of segmentation and transcript assignment in spatial omics data.
# https://github.com/LazDaria/SegTraQ

FROM python:3.12-slim

LABEL authors="Priyal Tripathi" \
    description="Docker image containing SegTraQ and its dependencies for segmentation QC"

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV NUMBA_CACHE_DIR='tmp'
ENV MPLCONFIGDIR='tmp/matplotlib'
ENV XDG_CACHE_HOME='tmp'

# Install system dependencies required by geopandas, rasterio, rtree, igraph
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgeos-dev \
    libgdal-dev \
    libspatialindex-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Install segtraq with pinned version and its key dependencies
RUN pip install --no-cache-dir \
    segtraq==0.0.3 \
    spatialdata>=0.7.2 \
    spatialdata-io>=0.1.4 \
    anndata>=0.12 \
    scanpy \
    squidpy>=1.6.2 \
    geopandas \
    igraph \
    rtree \
    rasterio \
    ovrlpy>=1.1.0
