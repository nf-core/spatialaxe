#!/usr/bin/env python3
"""
Combined Image QC Script

This script combines three image QC scripts into one:
- image_qc_roi_processing.py: GPU-accelerated tile/pixel-level focus maps
- image_qc_processing.py: Cell-based QC with Quarto figures
- image_qc_mapping_to_cells.py: Maps tile results to cells

Author: Hanneke Okkenhaug, Malwina Prater
"""

from __future__ import annotations

import gc
import logging
import math
import os
import sys
import time
import traceback
import click
import json
import numpy as np
import pandas as pd
import tifffile
import zarr
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
from napari_skimage_regionprops import regionprops_table
from skimage.segmentation import clear_border
from skimage import measure, color, morphology
import napari_simpleitk_image_processing as nsitk
import seaborn as sns
from sklearn.preprocessing import RobustScaler
from tifffile import imread
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any
from numpy.typing import NDArray
from scipy.ndimage import gaussian_laplace as scipy_gaussian_laplace
from scipy.ndimage import laplace as scipy_laplace
from scipy.ndimage import uniform_filter as scipy_uniform_filter
from scipy.ndimage import mean as ndimage_mean
from scipy import ndimage
from sklearn.mixture import GaussianMixture

import snr_metrics

# Set matplotlib to use a non-interactive backend
import matplotlib

matplotlib.use("Agg")

# GPU backend detection (CuPy)
try:
    import cupy as cp  # type: ignore[import-untyped]
    import cupyx.scipy.ndimage  # type: ignore[import-untyped]
    from cupyx.scipy.ndimage import laplace as cupy_laplace  # type: ignore[import-untyped]
    from cupyx.scipy.ndimage import uniform_filter as cupy_uniform_filter  # type: ignore[import-untyped]

    def cupy_gaussian_laplace(image, sigma):  # type: ignore[misc]
        """CuPy Laplacian of Gaussian: Gaussian smooth then Laplacian."""
        from cupyx.scipy.ndimage import gaussian_filter as _gf  # type: ignore[import-untyped]

        return cupy_laplace(_gf(image, sigma=sigma))

    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

# Xenium pixel size in micrometers (used for coordinate conversions)
XENIUM_PIXEL_SIZE_UM = 0.2125

# Default CCFS threshold for classifying cells as low nuclear texture quality.
# CCFS measures per-cell nuclear contrast (local_var/local_mean), NOT optical blur.
# Calibrated on 5 samples (2026-04-02): cells below 0.02 lose >50% transcripts
# relative to top-50% CCFS cells in lung/liver.  Brain tissue is confounded by
# nucleus size (large neurons score lower) — interpret with caution.
DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD = 0.02

#       sensitivity for blurred cell detection (catches more background regions)
ROI_INTENSITY_THRESHOLD = 100.0
ROI_MIN_TISSUE_COVERAGE_FOR_INTENSITY_QC = 0.5


# Fallback intensity critical thresholds (used when YAML not provided)
_INTENSITY_CRITICAL_DEFAULTS = {"dapi": 500, "boundary": 100, "intrna": 300}


def _load_qc_thresholds(yaml_path: str | None) -> dict:
    """Load thresholds from roi_image_qc_thresholds.yaml.

    Returns the ``image_qc`` sub-dict, or an empty dict if the file
    is missing / unreadable.
    """
    if yaml_path is None:
        return {}
    try:
        import yaml

        p = Path(yaml_path)
        if not p.exists():
            logging.warning("Thresholds YAML not found: %s", yaml_path)
            return {}
        with open(p, encoding="utf-8") as f:
            d = yaml.safe_load(f) or {}
        return d.get("image_qc", {}) if isinstance(d, dict) else {}
    except Exception as exc:
        logging.warning("Failed to load thresholds YAML: %s", exc)
        return {}


# Focus score percentile: Percentile of raw focus scores to use as threshold
# Applied to: tiles only after exclusion of low-intensity tiles (intensity >= ROI_INTENSITY_THRESHOLD)
# Purpose: Separate blurred vs in-focus tiles

ROI_FOCUS_SCORE_PERCENTILE = 5.0


def _run_figure_task(fn):
    """Wrapper for multiprocessing figure tasks that logs exceptions before exit.

    When using fork-based multiprocessing, child process tracebacks are lost
    if the child crashes. This wrapper catches and logs the full traceback
    in the child process before re-raising, so errors are visible in logs.
    """
    try:
        fn()
    except Exception:
        logging.error(f"Figure task {fn.__name__} failed:\n{traceback.format_exc()}")
        raise


def open_zarr(path: Path, zarr3: bool = False) -> zarr.Group:
    if zarr3:
        store = (
            zarr.storage.ZipStore(path)
            if path.suffix == ".zip"
            else zarr.storage.LocalStore(path)
        )
        return zarr.open_group(store=store, mode="r")
    else:
        """Open a Zarr file (compatible with zarr < 3)"""
        store = (
            zarr.ZipStore(path, mode="r")
            if path.suffix == ".zip"
            else zarr.DirectoryStore(path)
        )
        return zarr.group(store=store)


def load_and_prepare_data(xenium_bundle_dir, outdir):
    """Load all required data files and prepare paths"""

    xenium_bundle_dir = Path(xenium_bundle_dir)
    # Resolve outdir to absolute path to avoid nested directory issues
    # If outdir is already absolute, resolve() returns it as-is
    # If outdir is relative, resolve() makes it absolute relative to current working directory
    outdir = Path(outdir).resolve()

    # Required paths
    cells_parquet_path = xenium_bundle_dir / "cells.parquet"
    clusters_csv_path = (
        xenium_bundle_dir
        / "analysis"
        / "clustering"
        / "gene_expression_kmeans_10_clusters"
        / "clusters.csv"
    )
    umap_path = (
        xenium_bundle_dir
        / "analysis"
        / "umap"
        / "gene_expression_2_components"
        / "projection.csv"
    )
    cell_masks_path = xenium_bundle_dir / "cells.zarr.zip"
    morphology_focus_dir = xenium_bundle_dir / "morphology_focus"

    # Create output directories
    figures_dir = outdir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    return {
        "xenium_bundle_dir": xenium_bundle_dir,
        "outdir": outdir,
        "figures_dir": figures_dir,
        "cells_parquet_path": cells_parquet_path,
        "clusters_csv_path": clusters_csv_path,
        "umap_path": umap_path,
        "cell_masks_path": cell_masks_path,
        "morphology_focus_dir": morphology_focus_dir,
    }


def load_morphology_images(xoa_morphology_files):
    """Load downsampled (level 3) morphology images adaptively.

    Channel 0 (DAPI) is always loaded. Channels 1 (Boundary) and 2
    (Interior) are loaded only when the corresponding file is present,
    supporting both DAPI-only bundles (1 channel) and full
    DAPI/Boundary/Interior bundles (3 channels).

    Returns ``(small0, small1, small2)`` where ``small1`` and/or ``small2``
    may be ``None`` when the matching morphology file is absent.
    """

    small0 = tifffile.imread(
        xoa_morphology_files[0], is_ome=False, level=3, aszarr=False
    )

    small1 = None
    if len(xoa_morphology_files) > 1 and Path(xoa_morphology_files[1]).exists():
        small1 = tifffile.imread(
            xoa_morphology_files[1], is_ome=False, level=3, aszarr=False
        )

    small2 = None
    if len(xoa_morphology_files) > 2 and Path(xoa_morphology_files[2]).exists():
        small2 = tifffile.imread(
            xoa_morphology_files[2], is_ome=False, level=3, aszarr=False
        )

    return small0, small1, small2


def _load_morphology_channels(
    xoa_morphology_files,
    *,
    level: int = 0,
):
    """
    Load DAPI/Boundary/IntRNA channels robustly from either:
    - a single multi-channel OME-TIFF, or
    - separate single-channel files (e.g. *_0000, *_0001, *_0002).
    """

    def _split_channels(arr):
        # Support both channel-first (C,Y,X) and channel-last (Y,X,C).
        # .copy() on each slice so the parent 3D array can be GC'd.
        if arr.ndim == 2:
            return arr, None, None
        if arr.ndim != 3:
            raise ValueError(f"Unexpected morphology array shape: {arr.shape}")
        if arr.shape[0] <= 4 and arr.shape[1] > 16 and arr.shape[2] > 16:
            d = arr[0].copy()
            b = arr[1].copy() if arr.shape[0] > 1 else None
            r = arr[2].copy() if arr.shape[0] > 2 else None
            return d, b, r
        if arr.shape[2] <= 4 and arr.shape[0] > 16 and arr.shape[1] > 16:
            d = arr[:, :, 0].copy()
            b = arr[:, :, 1].copy() if arr.shape[2] > 1 else None
            r = arr[:, :, 2].copy() if arr.shape[2] > 2 else None
            return d, b, r
        # Conservative fallback: prefer first axis as channels.
        d = arr[0].copy()
        b = arr[1].copy() if arr.shape[0] > 1 else None
        r = arr[2].copy() if arr.shape[0] > 2 else None
        return d, b, r

    primary = tifffile.imread(
        xoa_morphology_files[0], is_ome=False, level=level, aszarr=False
    )

    if primary.ndim == 3:
        dapi, boundary, intrna = _split_channels(primary)
        return dapi, boundary, intrna

    dapi = primary
    boundary = None
    intrna = None

    if len(xoa_morphology_files) > 1 and Path(xoa_morphology_files[1]).exists():
        try:
            b = tifffile.imread(
                xoa_morphology_files[1], is_ome=False, level=level, aszarr=False
            )
            boundary = _split_channels(b)[0] if getattr(b, "ndim", 0) == 3 else b
        except Exception as e:
            logging.warning(
                "Boundary channel load failed for %s (continuing with DAPI): %s",
                xoa_morphology_files[1],
                e,
            )
            boundary = None
    if len(xoa_morphology_files) > 2 and Path(xoa_morphology_files[2]).exists():
        try:
            r = tifffile.imread(
                xoa_morphology_files[2], is_ome=False, level=level, aszarr=False
            )
            intrna = _split_channels(r)[0] if getattr(r, "ndim", 0) == 3 else r
        except Exception as e:
            logging.warning(
                "IntRNA channel load failed for %s (continuing with DAPI): %s",
                xoa_morphology_files[2],
                e,
            )
            intrna = None

    return dapi, boundary, intrna


class _LazyTiffChannel:
    """Lazy 2-D view into one channel page of a TIFF file.

    Supports ``[y0:y1, x0:x1]`` slicing.  For tiled TIFFs (production
    Xenium images) only the overlapping tiles are decoded, keeping I/O
    minimal.  For non-tiled TIFFs (e.g. test images written by
    ``tifffile.imwrite``) the full page is cached on first access.

    Thread-safe: a lock serialises file-handle reads so that
    ``_process_tile_on_gpu`` can call ``[slice]`` from multiple threads.
    """

    def __init__(self, page):
        """Accept a ``TiffPage`` or ``TiffFrame``.

        ``TiffFrame`` (non-first pages in a multi-page TIFF) lacks
        ``imagelength`` / ``is_tiled`` — we fall back to ``.shape`` and
        always use the cached-full-read path for frames.
        """
        import threading

        self._page = page
        # .shape works on both TiffPage and TiffFrame
        self.shape: tuple[int, int] = (page.shape[0], page.shape[1])
        self._lock = threading.Lock()
        self._cached_data: NDArray | None = None
        # TiffFrame has no is_tiled; treat as non-tiled (fallback path)
        self._is_tiled: bool = getattr(page, "is_tiled", False)

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def __getitem__(self, key):
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("_LazyTiffChannel supports only [y0:y1, x0:x1] slicing")
        yslice, xslice = key
        y0, y1, _ = yslice.indices(self.shape[0])
        x0, x1, _ = xslice.indices(self.shape[1])
        height = y1 - y0
        width = x1 - x0
        if height <= 0 or width <= 0:
            return np.empty((0, 0), dtype=self._page.dtype)
        if self._is_tiled:
            return self._read_region_tiled(y0, x0, height, width)
        return self._read_region_fallback(y0, x0, height, width)

    # ------------------------------------------------------------------
    # tiled path — decode only the tiles that overlap the request
    # ------------------------------------------------------------------

    def _read_region_tiled(self, y0: int, x0: int, height: int, width: int) -> NDArray:
        page = self._page
        tw, th = page.tilewidth, page.tilelength
        im_w, im_h = page.imagewidth, page.imagelength

        y1 = min(y0 + height, im_h)
        x1 = min(x0 + width, im_w)

        tile_y0 = y0 // th
        tile_x0 = x0 // tw
        tile_y1 = int(np.ceil(y1 / th))
        tile_x1 = int(np.ceil(x1 / tw))
        tiles_per_row = int(np.ceil(im_w / tw))

        buf_h = (tile_y1 - tile_y0) * th
        buf_w = (tile_x1 - tile_x0) * tw
        out = np.empty((buf_h, buf_w), dtype=page.dtype)

        fh = page.parent.filehandle
        with self._lock:
            for ti in range(tile_y0, tile_y1):
                for tj in range(tile_x0, tile_x1):
                    index = ti * tiles_per_row + tj
                    offset = page.dataoffsets[index]
                    bytecount = page.databytecounts[index]
                    fh.seek(offset)
                    data = fh.read(bytecount)
                    tile_arr, _indices, _shape = page.decode(
                        data, index, jpegtables=page.jpegtables
                    )
                    tile_2d = tile_arr.squeeze()
                    oi = (ti - tile_y0) * th
                    oj = (tj - tile_x0) * tw
                    out[oi : oi + th, oj : oj + tw] = tile_2d

        ry0 = y0 - tile_y0 * th
        rx0 = x0 - tile_x0 * tw
        return out[ry0 : ry0 + (y1 - y0), rx0 : rx0 + (x1 - x0)]

    # ------------------------------------------------------------------
    # fallback path — cache full page (non-tiled / test images)
    # ------------------------------------------------------------------

    def _read_region_fallback(
        self, y0: int, x0: int, height: int, width: int
    ) -> NDArray:
        with self._lock:
            if self._cached_data is None:
                self._cached_data = self._page.asarray()
        return self._cached_data[y0 : y0 + height, x0 : x0 + width]


def _open_morphology_lazy(
    xoa_morphology_files,
    *,
    level: int = 0,
) -> tuple[list, tuple[int, int]]:
    """Open morphology channels as lazy TIFF page wrappers (no pixel data loaded).

    Uses ``tifffile.TiffFile`` directly instead of the zarr store interface,
    avoiding the zarr v3 dependency introduced by ``tifffile.imread(aszarr=True)``.

    Args:
        xoa_morphology_files: List of OME-TIFF file paths.
        level: Resolution level to open (0 = full resolution).

    Returns:
        (channels, image_shape) where channels is [dapi, boundary, intrna]
        (None for missing channels) and image_shape is (height, width).
    """
    tiff_handles: list[tifffile.TiffFile] = []  # prevent GC

    primary_tif = tifffile.TiffFile(str(xoa_morphology_files[0]))
    tiff_handles.append(primary_tif)

    # Use PHYSICAL pages in the file, not OME series pages.
    # Xenium multi-file OME-TIFFs report series.shape=(4,H,W) referencing
    # all files, but each file has only 1 physical page. series.pages would
    # include stubs for data in other files that decode to zeros.
    n_physical = len(primary_tif.pages)
    if level > 0:
        # For pyramid levels, use the series API to navigate levels
        series = primary_tif.series[0]
        if level < len(series.levels):
            pages = list(series.levels[level].pages)
        else:
            pages = list(primary_tif.pages)
    else:
        pages = list(primary_tif.pages)

    n_pages = len(pages)

    # Multi-channel single file: truly multi-page TIFF (each page = channel)
    # Only enter this path if the file physically contains multiple pages
    if n_pages >= 2 and n_physical >= 2:
        # Each page is one channel (C, H, W layout across pages)
        dapi = _LazyTiffChannel(pages[0])
        shape = dapi.shape
        boundary = _LazyTiffChannel(pages[1]) if n_pages > 1 else None
        intrna = _LazyTiffChannel(pages[2]) if n_pages > 2 else None
        # Attach TiffFile handles to prevent garbage collection
        dapi._tiff_handles = tiff_handles  # type: ignore[attr-defined]
        return [dapi, boundary, intrna], shape

    # Single-channel (or single-page) primary
    # Check if the single page is actually multi-channel (C, H, W) in one page
    page0 = pages[0]
    arr_shape = page0.shape  # could be (H, W) or (C, H, W)
    if len(arr_shape) == 3:
        # Multi-channel packed into one page — fall back to cached full read
        full = page0.asarray()
        if arr_shape[0] <= 4 and arr_shape[1] > 16 and arr_shape[2] > 16:
            # Channel-first (C, H, W)
            shape = (arr_shape[1], arr_shape[2])
            ch_axis = 0
        elif arr_shape[2] <= 4 and arr_shape[0] > 16 and arr_shape[1] > 16:
            # Channel-last (H, W, C)
            shape = (arr_shape[0], arr_shape[1])
            ch_axis = 2
        else:
            shape = (arr_shape[1], arr_shape[2])
            ch_axis = 0
        n_ch = arr_shape[ch_axis]
        dapi = _LazyTiffChannel.__new__(_LazyTiffChannel)
        dapi._page = page0
        dapi.shape = shape
        import threading

        dapi._lock = threading.Lock()
        dapi._cached_data = np.take(full, 0, axis=ch_axis)
        boundary_ch = _LazyTiffChannel.__new__(_LazyTiffChannel) if n_ch > 1 else None
        if boundary_ch is not None:
            boundary_ch._page = page0
            boundary_ch.shape = shape
            boundary_ch._lock = threading.Lock()
            boundary_ch._cached_data = np.take(full, 1, axis=ch_axis)
        intrna_ch = _LazyTiffChannel.__new__(_LazyTiffChannel) if n_ch > 2 else None
        if intrna_ch is not None:
            intrna_ch._page = page0
            intrna_ch.shape = shape
            intrna_ch._lock = threading.Lock()
            intrna_ch._cached_data = np.take(full, 2, axis=ch_axis)
        dapi._tiff_handles = tiff_handles  # type: ignore[attr-defined]
        return [dapi, boundary_ch, intrna_ch], shape

    # Truly single-channel 2-D page
    dapi = _LazyTiffChannel(page0)
    shape = dapi.shape
    boundary = None
    intrna = None

    if len(xoa_morphology_files) > 1 and Path(xoa_morphology_files[1]).exists():
        try:
            t = tifffile.TiffFile(str(xoa_morphology_files[1]))
            tiff_handles.append(t)
            boundary = _LazyTiffChannel(t.pages[0])
        except Exception as e:
            logging.warning(
                "Boundary TIFF open failed for %s: %s", xoa_morphology_files[1], e
            )

    if len(xoa_morphology_files) > 2 and Path(xoa_morphology_files[2]).exists():
        try:
            t = tifffile.TiffFile(str(xoa_morphology_files[2]))
            tiff_handles.append(t)
            intrna = _LazyTiffChannel(t.pages[0])
        except Exception as e:
            logging.warning(
                "IntRNA TIFF open failed for %s: %s", xoa_morphology_files[2], e
            )

    # Attach TiffFile handles to prevent garbage collection
    dapi._tiff_handles = tiff_handles  # type: ignore[attr-defined]
    return [dapi, boundary, intrna], shape


def generate_tissue_mask(
    xoa_morphology_files,
    small0,
    small1,
    small2,
    threshold_percentile=60,
    min_size_edge=500000,
    min_size_hole=1500,
    dense_intensity_region_percentile=97,
    min_size_dense_intensity_region=500,
    downsample_factor=8,
):
    """
    Generate tissue masks and distance maps from morphology images (cell/segmentation independent).

    Parameters:
    -----------
    xoa_morphology_files : list
        List of paths to morphology image files (for compatibility, not used if small0/1/2 provided)
    small0, small1, small2 : numpy.ndarray
        Downsampled morphology images (DAPI, Boundary, Interior)
    threshold_percentile : float, optional
        Percentile for thresholding (default: 60)
    min_size_edge : int, optional
        Minimum size for edge objects in downsampled space (default: 500000)
    min_size_hole : int, optional
        Minimum size for holes in downsampled space (default: 1500)
    dense_intensity_region_percentile : float, optional
        Percentile for dense intensity region detection (default: 97)
    min_size_dense_intensity_region : int, optional
        Minimum size for dense intensity regions (default: 500)
    downsample_factor : int, optional
        Downsampling factor (default: 8, for level 3)

    Returns:
    --------
    tuple
        (whole_sample, holes, dense_intensity_regions, distance_map, distance_map2)
        - whole_sample: Labeled tissue mask
        - holes: Labeled holes mask
        - dense_intensity_regions: Labeled dense intensity regions mask
        - distance_map: Distance to edge map
        - distance_map2: Distance to nearest hole map
    """
    # Intensity-based thresholding
    t = np.percentile(small0, threshold_percentile)
    thresh1 = small0 < t
    thresh2 = small0 >= t

    # Find the edge of the sample
    objects = measure.label(thresh1)
    mask = morphology.remove_small_objects(objects, min_size=min_size_edge)
    edge_sample = measure.label(mask)
    distance_map = nsitk.signed_maurer_distance_map(edge_sample)

    # Find the holes in the sample
    noborder = clear_border(objects)
    holes = morphology.remove_small_objects(noborder, min_size=min_size_hole)
    small_objects = noborder ^ holes
    test_mask = measure.label(thresh2) + small_objects
    sample = test_mask > 1
    whole_sample = measure.label(sample)
    distance_map2 = nsitk.signed_maurer_distance_map(holes)

    # Detecting dense intensity regions
    # Boundary (small1) / Interior (small2) channels are optional: DAPI-only
    # bundles pass them as None, in which case only the DAPI channel
    # contributes to the dense-intensity-region detection.
    t0 = np.percentile(small0, dense_intensity_region_percentile)
    thresh_sum = small0 >= t0
    if small1 is not None:
        t1 = np.percentile(small1, dense_intensity_region_percentile)
        thresh_sum = thresh_sum + (small1 >= t1)
    if small2 is not None:
        t2 = np.percentile(small2, dense_intensity_region_percentile)
        thresh_sum = thresh_sum + (small2 >= t2)
    thresh_fill = nsitk.binary_fill_holes(thresh_sum)
    objects_art = measure.label(thresh_fill)
    dense_intensity_regions = morphology.remove_small_objects(
        objects_art, min_size=min_size_dense_intensity_region
    )

    return whole_sample, holes, dense_intensity_regions, distance_map, distance_map2


# ---------------------------------------------------------------------------
# Internal helpers (from focus_score_maps)
# ---------------------------------------------------------------------------

_EPSILON: float = 1e-8


def _get_backend(
    use_gpu: bool,
) -> tuple[Any, Any, Any]:
    """Return the array module, uniform_filter, and laplace for the chosen backend.

    Args:
        use_gpu: If ``True`` and CuPy is available, return GPU primitives.

    Returns:
        Tuple of ``(array_module, uniform_filter_fn, laplace_fn)``.

    Raises:
        RuntimeError: If ``use_gpu`` is ``True`` but CuPy is not installed.
    """
    if use_gpu:
        if not HAS_CUPY:
            raise RuntimeError(
                "use_gpu=True but CuPy is not installed. "
                "Install CuPy or set use_gpu=False."
            )
        return cp, cupy_uniform_filter, cupy_laplace
    return np, scipy_uniform_filter, scipy_laplace


def _to_device(image: NDArray[np.generic], xp: Any) -> Any:
    """Move a numpy array to the target device (no-op for NumPy).

    Args:
        image: Input numpy array.
        xp: Array module (``numpy`` or ``cupy``).

    Returns:
        Array on the target device.
    """
    if xp is np:
        return image
    return cp.asarray(image)


def _to_numpy(arr: Any, xp: Any) -> NDArray[np.float32]:
    """Move an array back to host memory as float32 numpy (no-op for NumPy).

    Args:
        arr: Array on device.
        xp: Array module that created *arr*.

    Returns:
        numpy float32 array on the host.
    """
    if xp is np:
        return np.asarray(arr, dtype=np.float32)
    return cp.asnumpy(arr).astype(np.float32)


def _sanitize(arr: NDArray[np.float32]) -> NDArray[np.float32]:
    """Replace NaN and Inf values with zero.

    Args:
        arr: Input array (modified in-place).

    Returns:
        The same array with non-finite values set to zero.
    """
    arr[~np.isfinite(arr)] = 0.0
    return arr


def detect_gpu_ids() -> list[int]:
    """Detect available CUDA GPUs.

    Returns:
        List of GPU device IDs. Empty list if CuPy is not available or no GPUs
        are detected.
    """
    if not HAS_CUPY:
        return []
    try:
        n_devices = cp.cuda.runtime.getDeviceCount()
        return list(range(n_devices))
    except Exception:
        return []


def _compute_adaptive_tile_size(
    height: int,
    width: int,
    n_gpus: int,
    gpu_mem_bytes: int | None = None,
    target_utilization: float = 0.65,
    min_tiles_per_gpu: int = 2,
) -> int:
    """Compute tile size that maximizes GPU memory utilization.

    Balances two constraints:
    1. Each tile's peak GPU memory stays within target_utilization of VRAM.
    2. Total tiles >= n_gpus * min_tiles_per_gpu for load balancing.

    Args:
        height: Image height in pixels.
        width: Image width in pixels.
        n_gpus: Number of available GPUs.
        gpu_mem_bytes: Total GPU VRAM in bytes. Auto-detected if None.
        target_utilization: Fraction of VRAM to target per tile (default 0.65).
        min_tiles_per_gpu: Minimum tiles per GPU for load balancing.

    Returns:
        Tile size in pixels (square tiles).
    """
    if gpu_mem_bytes is None:
        gpu_mem_bytes = cp.cuda.Device(0).mem_info[1]

    # Peak memory per pixel: ~6 concurrent float32 arrays
    BYTES_PER_PIXEL = 6 * 4  # 24 bytes

    # Max tile dim from GPU memory constraint
    max_pixels = int(gpu_mem_bytes * target_utilization / BYTES_PER_PIXEL)
    max_tile_dim = int(math.sqrt(max_pixels))

    # Min tiles needed for load balancing
    min_tiles = max(n_gpus * min_tiles_per_gpu, 1)

    # Minimum useful tile dimension (must fit the convolution window)
    min_dim = 256

    # Compute tile_size that produces at least min_tiles
    # Start from max_tile_dim and shrink until we have enough tiles
    tile_size = max_tile_dim
    while tile_size > min_dim:
        n_y = math.ceil(height / tile_size)
        n_x = math.ceil(width / tile_size)
        if n_y * n_x >= min_tiles:
            break
        tile_size = int(tile_size * 0.7)  # shrink by 30%

    # Clamp to minimum useful dimension
    tile_size = max(min_dim, tile_size)

    # If image is smaller than tile_size, just use image dims
    if height <= tile_size and width <= tile_size:
        tile_size = max(height, width)

    return tile_size


def _compute_tile_grid(
    height: int,
    width: int,
    tile_size: int = 8192,
    overlap: int = 17,
) -> list[dict[str, int]]:
    """Compute overlapping tile coordinates for tiled convolution.

    Each tile has a "write" region (non-overlapping, covers full image) and
    a "read" region (expanded by overlap, clipped to image bounds).

    Args:
        height: Image height in pixels.
        width: Image width in pixels.
        tile_size: Core tile dimension (before overlap).
        overlap: Border pixels to add for convolution safety.

    Returns:
        List of tile spec dicts with keys: read_y0, read_y1, read_x0, read_x1,
        write_y0, write_y1, write_x0, write_x1, trim_top, trim_bottom,
        trim_left, trim_right.
    """
    tiles = []
    for y in range(0, height, tile_size):
        for x in range(0, width, tile_size):
            wy0 = y
            wy1 = min(y + tile_size, height)
            wx0 = x
            wx1 = min(x + tile_size, width)

            ry0 = max(0, wy0 - overlap)
            ry1 = min(height, wy1 + overlap)
            rx0 = max(0, wx0 - overlap)
            rx1 = min(width, wx1 + overlap)

            tiles.append(
                {
                    "read_y0": ry0,
                    "read_y1": ry1,
                    "read_x0": rx0,
                    "read_x1": rx1,
                    "write_y0": wy0,
                    "write_y1": wy1,
                    "write_x0": wx0,
                    "write_x1": wx1,
                    "trim_top": wy0 - ry0,
                    "trim_bottom": ry1 - wy1,
                    "trim_left": wx0 - rx0,
                    "trim_right": rx1 - wx1,
                }
            )
    return tiles


# ---------------------------------------------------------------------------
# Focus-map computation
# ---------------------------------------------------------------------------


def compute_ccfs_map(
    image: NDArray[np.generic],
    window_size: int = 35,
    use_gpu: bool = False,
    gpu_id: int = 0,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Compute a per-pixel CCFS (Coefficient of Contrast Focus Score) map.

    The CCFS focus score at each pixel is defined as::

        focus = local_var / local_mean

    which is equivalent to ``std**2 / mean`` computed over a square window
    centred on that pixel.

    Args:
        image: 2-D input image (any numeric dtype).
        window_size: Side length of the square averaging window.
        use_gpu: Use CuPy GPU backend when ``True``.
        gpu_id: CUDA device ID to use when ``use_gpu=True``.

    Returns:
        Tuple of ``(focus_map, mean_map)`` — both float32 arrays with the
        same shape as *image*.

    Raises:
        ValueError: If *image* is not 2-D or is smaller than the window in
            either dimension.
        RuntimeError: If *use_gpu* is ``True`` but CuPy is unavailable.
    """
    if image.ndim != 2:
        raise ValueError(f"Expected a 2-D image, got shape {image.shape}.")
    if image.shape[0] < window_size or image.shape[1] < window_size:
        raise ValueError(
            f"Image shape {image.shape} is smaller than window_size "
            f"{window_size} in at least one dimension."
        )

    xp, uniform_filter, _ = _get_backend(use_gpu)

    if use_gpu:
        with cp.cuda.Device(gpu_id):
            image_f = cp.asarray(image.astype(np.float32, copy=False))
            local_mean = uniform_filter(image_f, size=window_size)
            local_sq_mean = uniform_filter(image_f**2, size=window_size)
            # Promote to float64 for subtraction to avoid catastrophic cancellation
            local_var = (
                local_sq_mean.astype(cp.float64) - local_mean.astype(cp.float64) ** 2
            )
            local_var = xp.maximum(local_var, 0.0)
            focus_map_dev = (
                local_var / (local_mean.astype(cp.float64) + _EPSILON)
            ).astype(cp.float32)
            focus_map = _to_numpy(focus_map_dev, xp)
            mean_map = _to_numpy(local_mean, xp)
    else:
        image_f = image.astype(np.float32, copy=False)
        local_mean = uniform_filter(image_f, size=window_size)
        sq = image_f**2
        del image_f
        local_sq_mean = uniform_filter(sq, size=window_size)
        del sq
        # Promote to float64 for the subtraction to avoid catastrophic
        # cancellation on high-intensity uint16 images.
        # Explicit steps to limit peak memory (avoid 3+ simultaneous f64 temps).
        local_var = local_sq_mean.astype(np.float64)
        del local_sq_mean
        temp_mean_f64 = local_mean.astype(np.float64)
        local_var -= temp_mean_f64**2
        del temp_mean_f64
        local_var = np.maximum(local_var, 0.0)
        # Compute focus_map = local_var / (local_mean + eps).
        # Save mean_map first, then free local_mean before creating f64 temp.
        mean_map = local_mean.astype(np.float32)
        temp_denom = local_mean.astype(np.float64)
        del local_mean
        temp_denom += _EPSILON
        local_var /= temp_denom
        del temp_denom
        focus_map = local_var.astype(np.float32)
        del local_var

    return _sanitize(focus_map), _sanitize(mean_map)


def compute_laplacian_variance_map(
    image: NDArray[np.generic],
    window_size: int = 35,
    use_gpu: bool = False,
    gpu_id: int = 0,
    lap_sigma: float = 1.0,
) -> NDArray[np.float32]:
    """Compute a per-pixel windowed Laplacian-variance focus map.

    This is a *local* variant of the standard Laplacian variance focus metric
    (Pech-Pacheco et al., 2000).  Instead of computing a single global
    variance, it produces a spatial map where each pixel holds the variance
    of the Laplacian of Gaussian (LoG) response inside the surrounding
    *window_size* × *window_size* neighbourhood::

        lap_var = uniform_filter(lap**2) - uniform_filter(lap)**2

    Using LoG (``gaussian_laplace`` with *lap_sigma*) rather than the bare
    3×3 Laplacian suppresses pixel-level noise and makes the metric more
    specific to genuine edge content (Sun et al., 2004; Pertuz et al., 2013).

    Variance is computed in float64 to avoid catastrophic cancellation in the
    ``E[X²] − E[X]²`` formula, then cast back to float32.

    Args:
        image: 2-D input image (any numeric dtype).
        window_size: Side length of the square averaging window.
        use_gpu: Use CuPy GPU backend when ``True``.
        gpu_id: CUDA device ID to use when ``use_gpu=True``.
        lap_sigma: Gaussian sigma for LoG pre-smoothing (pixels).
            Set to 0 to revert to the bare Laplacian.

    Returns:
        Float32 array with the same shape as *image*.

    Raises:
        ValueError: If *image* is not 2-D or is smaller than the window in
            either dimension.
        RuntimeError: If *use_gpu* is ``True`` but CuPy is unavailable.
    """
    if image.ndim != 2:
        raise ValueError(f"Expected a 2-D image, got shape {image.shape}.")
    if image.shape[0] < window_size or image.shape[1] < window_size:
        raise ValueError(
            f"Image shape {image.shape} is smaller than window_size "
            f"{window_size} in at least one dimension."
        )

    xp, uniform_filter, laplace_fn = _get_backend(use_gpu)

    if use_gpu:
        with cp.cuda.Device(gpu_id):
            image_f = cp.asarray(image.astype(np.float32, copy=False))
            if lap_sigma > 0:
                lap = cupy_gaussian_laplace(image_f, sigma=lap_sigma)
            else:
                lap = laplace_fn(image_f)
            lap = lap.astype(cp.float64)
            lap_mean = uniform_filter(lap, size=window_size)
            lap_sq_mean = uniform_filter(lap * lap, size=window_size)
            lap_var = lap_sq_mean - lap_mean**2
            lap_var = xp.maximum(lap_var, 0.0)
            result = _to_numpy(lap_var.astype(cp.float32), xp)
    else:
        image_f = image.astype(np.float32, copy=False)
        if lap_sigma > 0:
            lap = scipy_gaussian_laplace(image_f, sigma=lap_sigma).astype(np.float64)
        else:
            lap = laplace_fn(image_f).astype(np.float64)
        del image_f
        lap_sq = lap * lap
        lap_mean = uniform_filter(lap, size=window_size)
        del lap
        lap_sq_mean = uniform_filter(lap_sq, size=window_size)
        del lap_sq
        np.square(lap_mean, out=lap_mean)  # in-place to avoid temporary
        lap_var = lap_sq_mean - lap_mean
        del lap_sq_mean, lap_mean
        lap_var = np.maximum(lap_var, 0.0)
        result = lap_var.astype(np.float32)
        del lap_var

    return _sanitize(result)


def _compute_channel_maps_on_gpu(
    channel: NDArray[np.generic],
    window_size: int,
    gpu_id: int,
    include_laplacian: bool = False,
    lap_sigma: float = 1.0,
) -> dict[str, NDArray[np.float32]]:
    """Compute focus + mean maps for a single channel on a specific GPU.

    This fused implementation opens a single device context and transfers the
    image to the GPU once, computing CCFS (focus_map, mean_map) and optionally
    the Laplacian variance map within the same context. This avoids redundant
    CPU-to-GPU transfers of the full image.

    Args:
        channel: 2-D image array.
        window_size: Convolution window size.
        gpu_id: CUDA device ID.
        include_laplacian: Also compute Laplacian variance map.
        lap_sigma: Gaussian sigma for LoG pre-smoothing (0 = bare Laplacian).

    Returns:
        Dict with ``focus_map``, ``mean_map``, and optionally ``lap_var_map``.
    """
    if not HAS_CUPY:
        raise RuntimeError(
            "_compute_channel_maps_on_gpu requires CuPy but it is not installed."
        )

    pool = cp.get_default_memory_pool()

    with cp.cuda.Device(gpu_id):
        # Upload image to GPU
        image_f = cp.asarray(channel.astype(np.float32, copy=False))

        # --- CCFS (focus_map + mean_map) ---
        # Peak memory: max 3 arrays alive at once (~10.5GB for 34K×26K images)
        local_mean = cupy_uniform_filter(image_f, size=window_size)
        sq = image_f**2
        # Free image_f before allocating local_sq_mean to cap at 3 arrays
        del image_f
        pool.free_all_blocks()
        local_sq_mean = cupy_uniform_filter(sq, size=window_size)
        del sq
        # Promote to float64 for subtraction to avoid catastrophic cancellation
        local_var = (
            local_sq_mean.astype(cp.float64) - local_mean.astype(cp.float64) ** 2
        )
        del local_sq_mean
        cp.maximum(local_var, 0.0, out=local_var)
        focus_map_dev = (local_var / (local_mean.astype(cp.float64) + _EPSILON)).astype(
            cp.float32
        )
        del local_var

        # Transfer CCFS results to CPU immediately
        focus_map = cp.asnumpy(focus_map_dev).astype(np.float32)
        mean_map = cp.asnumpy(local_mean).astype(np.float32)
        del focus_map_dev, local_mean
        pool.free_all_blocks()

        # --- Laplacian variance (optional, re-upload image from CPU) ---
        result: dict[str, NDArray[np.float32]] = {
            "focus_map": _sanitize(focus_map),
            "mean_map": _sanitize(mean_map),
        }
        if include_laplacian:
            image_f = cp.asarray(channel.astype(np.float32, copy=False))
            if lap_sigma > 0:
                lap = cupy_gaussian_laplace(image_f, sigma=lap_sigma)
            else:
                lap = cupy_laplace(image_f)
            del image_f
            pool.free_all_blocks()
            lap = lap.astype(cp.float64)
            lap_mean = cupy_uniform_filter(lap, size=window_size)
            sq = lap * lap
            del lap
            pool.free_all_blocks()
            lap_sq_mean = cupy_uniform_filter(sq, size=window_size)
            del sq
            lap_var = lap_sq_mean - lap_mean**2
            del lap_sq_mean, lap_mean
            cp.maximum(lap_var, 0.0, out=lap_var)
            lap_var_np = cp.asnumpy(lap_var).astype(np.float32)
            del lap_var
            pool.free_all_blocks()
            result["lap_var_map"] = _sanitize(lap_var_np)

    return result


def _process_tile_on_gpu(
    channel_data,
    tile_spec: dict[str, int],
    window_size: int,
    gpu_id: int,
    include_laplacian: bool = False,
    lap_sigma: float = 1.0,
) -> dict[str, NDArray[np.float32]]:
    """Read a tile from an array-like, compute focus maps on GPU, trim overlap.

    Args:
        channel_data: Array-like supporting slicing (_LazyTiffChannel or numpy array).
        tile_spec: Dict from _compute_tile_grid() with read/write/trim keys.
        window_size: Convolution window size.
        gpu_id: CUDA device ID.
        include_laplacian: Also compute Laplacian variance map.
        lap_sigma: Gaussian sigma for LoG.

    Returns:
        Dict with trimmed 'focus_map', 'mean_map', and optionally 'lap_var_map'.
    """
    # Read tile (triggers actual I/O for lazy TIFF-backed arrays)
    tile = np.asarray(
        channel_data[
            tile_spec["read_y0"] : tile_spec["read_y1"],
            tile_spec["read_x0"] : tile_spec["read_x1"],
        ]
    )

    # Compute on GPU using existing function
    result = _compute_channel_maps_on_gpu(
        tile, window_size, gpu_id, include_laplacian, lap_sigma
    )

    # Trim overlap from each output array
    tt = tile_spec["trim_top"]
    tb = tile_spec["trim_bottom"]
    tl = tile_spec["trim_left"]
    tr = tile_spec["trim_right"]

    trimmed = {}
    for key, arr in result.items():
        h, w = arr.shape
        y1 = h - tb if tb > 0 else h
        x1 = w - tr if tr > 0 else w
        trimmed[key] = arr[tt:y1, tl:x1]
    return trimmed


def _compute_channel_maps_tiled(
    channel_data,
    image_shape: tuple[int, int],
    window_size: int,
    gpu_ids: list[int],
    include_laplacian: bool = False,
    lap_sigma: float = 1.0,
    gpu_mem_bytes: int | None = None,
) -> dict[str, "NDArray[np.float32] | None"]:
    """Compute focus maps for one channel using tiled multi-GPU processing.

    Tiles the image with convolution-safe overlap, distributes tiles across
    GPUs via ThreadPoolExecutor, and assembles results into full-res arrays.
    Tile size is computed adaptively to maximize GPU memory utilization.

    Args:
        channel_data: Array-like (_LazyTiffChannel or numpy) with shape matching image_shape.
        image_shape: (height, width).
        window_size: Convolution window size.
        gpu_ids: List of CUDA device IDs.
        include_laplacian: Also compute Laplacian variance.
        lap_sigma: Gaussian sigma for LoG.
        gpu_mem_bytes: Total GPU VRAM in bytes. Auto-detected if None.

    Returns:
        Dict with 'focus_map', 'mean_map', and optionally 'lap_var_map'
        (all full-resolution float32 arrays).
    """
    H, W = image_shape
    n_gpus = len(gpu_ids)
    # Overlap must cover both uniform_filter radius (window_size // 2) and the
    # Gaussian pre-smoothing in the Laplacian path (~3 * lap_sigma).  Using the
    # full window_size is safe for all kernels and adds negligible I/O overhead.
    overlap = window_size

    # Adaptive tile sizing: maximize GPU memory utilization
    tile_size = _compute_adaptive_tile_size(H, W, n_gpus, gpu_mem_bytes)
    tiles = _compute_tile_grid(H, W, tile_size, overlap)

    logging.info(
        f"    Tiled processing: {len(tiles)} tiles ({tile_size}px), "
        f"{n_gpus} GPU(s), overlap={overlap}px"
    )

    # Ensure CUDA_PATH is set for CuPy NVRTC kernel compilation in worker
    # threads.  CuPy auto-detects CUDA in the main thread but worker threads
    # can fail with "Failed to auto-detect CUDA root directory".  The pip
    # package nvidia-cuda-runtime-cu12 installs headers under
    # site-packages/nvidia/cuda_runtime/include/.
    if "CUDA_PATH" not in os.environ:
        try:
            import nvidia.cuda_runtime as _cr

            _cr_dir = _cr.__path__[0]  # .../site-packages/nvidia/cuda_runtime
            if os.path.isdir(os.path.join(_cr_dir, "include")):
                os.environ["CUDA_PATH"] = _cr_dir
                logging.info(f"    Set CUDA_PATH={_cr_dir}")
        except (ImportError, IndexError, AttributeError):
            pass  # CUDA_PATH remains unset; CuPy will try its own detection

    # Warm up CuPy kernel cache in the main thread.  NVRTC compiles kernels on
    # first use; running a tiny operation here populates the disk cache so that
    # worker threads hit the cache instead of compiling in parallel.
    try:
        _warmup = cp.array([1.0], dtype=cp.float32)
        cupyx.scipy.ndimage.uniform_filter(_warmup.reshape(1, 1), size=1)
        del _warmup
    except Exception:
        pass

    # Pre-allocate output arrays
    focus_out = np.empty((H, W), dtype=np.float32)
    mean_out = np.empty((H, W), dtype=np.float32)
    lap_out = np.empty((H, W), dtype=np.float32) if include_laplacian else None

    # Pipelined I/O + GPU: each worker reads its tile from TIFF then computes
    # on its assigned GPU.  The TIFF file handle lock serializes reads, but I/O
    # for tile N+1 overlaps with GPU compute for tile N (pipelining).  This is
    # faster than pre-reading all tiles because it allows I/O-GPU overlap.
    # True multi-GPU parallelism would require ProcessPoolExecutor (separate
    # GIL per process) or CuPy async streams.
    with ThreadPoolExecutor(max_workers=n_gpus) as executor:
        futures = {}
        for i, tile_spec in enumerate(tiles):
            assigned_gpu = gpu_ids[i % n_gpus]
            future = executor.submit(
                _process_tile_on_gpu,
                channel_data,
                tile_spec,
                window_size,
                assigned_gpu,
                include_laplacian,
                lap_sigma,
            )
            futures[future] = tile_spec

        for future in as_completed(futures):
            tile_spec = futures[future]
            result = future.result()
            wy0, wy1 = tile_spec["write_y0"], tile_spec["write_y1"]
            wx0, wx1 = tile_spec["write_x0"], tile_spec["write_x1"]
            focus_out[wy0:wy1, wx0:wx1] = result["focus_map"]
            mean_out[wy0:wy1, wx0:wx1] = result["mean_map"]
            if lap_out is not None and "lap_var_map" in result:
                lap_out[wy0:wy1, wx0:wx1] = result["lap_var_map"]

    focus_maps: dict[str, "np.ndarray | None"] = {
        "focus_map": _sanitize(focus_out),
        "mean_map": _sanitize(mean_out),
        "lap_var_map": _sanitize(lap_out) if lap_out is not None else None,
    }
    return focus_maps


def compute_all_focus_maps(
    channels: NDArray[np.generic],
    window_size: int = 35,
    use_gpu: bool = False,
    gpu_ids: list[int] | None = None,
    lap_sigma: float = 1.0,
) -> dict[str, NDArray[np.float32] | None]:
    """Compute focus maps for all available image channels.

    When ``gpu_ids`` contains multiple IDs, channel computations are
    distributed across GPUs in parallel using a thread pool. Each channel's
    convolution work runs entirely on one GPU; with 3 channels and 3+ GPUs
    all channels are processed concurrently.

    Args:
        channels: Either a 2-D array (single DAPI channel, shape ``(H, W)``)
            or a 3-D array with channels first (shape ``(C, H, W)``).
        window_size: Side length of the square averaging window.
        use_gpu: Use CuPy GPU backend when ``True``. Ignored when *gpu_ids*
            is provided (GPU is assumed).
        gpu_ids: List of CUDA device IDs for multi-GPU parallelism. If
            ``None`` and ``use_gpu=True``, uses device 0 only. If ``None``
            and ``use_gpu=False``, runs on CPU.

    Returns:
        Dictionary with the following keys (values are ``None`` when the
        corresponding channel is not present in *channels*):

        - ``dapi_focus_map`` — CCFS focus map for DAPI (channel 0).
        - ``dapi_mean_map`` — Local mean map for DAPI.
        - ``dapi_lap_var_map`` — Laplacian variance map for DAPI.
        - ``boundary_focus_map`` — CCFS focus map for Boundary (channel 1).
        - ``boundary_mean_map`` — Local mean map for Boundary.
        - ``intrna_focus_map`` — CCFS focus map for IntRNA (channel 2).
        - ``intrna_mean_map`` — Local mean map for IntRNA.

    Raises:
        ValueError: If *channels* has fewer than 2 or more than 3 dimensions.
    """
    # Resolve GPU configuration
    if gpu_ids is not None and len(gpu_ids) > 0:
        use_gpu = True
    elif use_gpu and gpu_ids is None:
        gpu_ids = [0]

    # Split channels — .copy() so the parent 3-D array can be freed.
    if channels.ndim == 2:
        dapi = channels.copy()
        boundary = None
        intrna = None
    elif channels.ndim == 3:
        n_ch = channels.shape[0]
        dapi = channels[0].copy()
        boundary = channels[1].copy() if n_ch > 1 else None
        intrna = channels[2].copy() if n_ch > 2 else None
    else:
        raise ValueError(
            f"Expected 2-D or 3-D array, got {channels.ndim}-D "
            f"(shape {channels.shape})."
        )
    del channels

    # ---- Multi-GPU path: distribute channels across GPUs ----
    if use_gpu and gpu_ids is not None and len(gpu_ids) > 0:
        # Build work items: (channel_name, channel_data, include_laplacian)
        work_items: list[tuple[str, NDArray[np.generic], bool]] = [
            ("dapi", dapi, True),  # DAPI always gets Laplacian
        ]
        if boundary is not None:
            work_items.append(("boundary", boundary, False))
        if intrna is not None:
            work_items.append(("intrna", intrna, False))

        # Assign GPUs round-robin
        results: dict[str, dict[str, NDArray[np.float32]]] = {}
        n_gpus = len(gpu_ids)
        logging.info(
            f"    Distributing {len(work_items)} channel(s) across {n_gpus} GPU(s): {gpu_ids}"
        )

        t_gpu = time.perf_counter()
        gpu_failed = False
        try:
            with ThreadPoolExecutor(max_workers=len(work_items)) as executor:
                futures = {}
                for idx, (name, ch_data, inc_lap) in enumerate(work_items):
                    assigned_gpu = gpu_ids[idx % n_gpus]
                    future = executor.submit(
                        _compute_channel_maps_on_gpu,
                        ch_data,
                        window_size,
                        assigned_gpu,
                        inc_lap,
                        lap_sigma,
                    )
                    futures[future] = name

                for future in as_completed(futures):
                    ch_name = futures[future]
                    results[ch_name] = future.result()
            logging.info(
                f"  [TIMING] GPU compute (multi-GPU, {len(work_items)} channels): {time.perf_counter() - t_gpu:.1f}s"
            )
        except Exception as gpu_err:
            logging.warning(
                f"  WARNING: GPU computation failed ({gpu_err}), falling back to CPU..."
            )
            gpu_failed = True

        if not gpu_failed:
            # Assemble output dict from GPU results
            dapi_res = results["dapi"]
            return {
                "dapi_focus_map": dapi_res["focus_map"],
                "dapi_mean_map": dapi_res["mean_map"],
                "dapi_lap_var_map": dapi_res.get("lap_var_map"),
                "boundary_focus_map": results["boundary"]["focus_map"]
                if "boundary" in results
                else None,
                "boundary_mean_map": results["boundary"]["mean_map"]
                if "boundary" in results
                else None,
                "intrna_focus_map": results["intrna"]["focus_map"]
                if "intrna" in results
                else None,
                "intrna_mean_map": results["intrna"]["mean_map"]
                if "intrna" in results
                else None,
            }
        # else: fall through to CPU path below

    # ---- Single-GPU or CPU fallback path ----
    # Reached when: no multi-GPU available, use_gpu=False, or GPU OOM fallback
    use_gpu = (
        False  # Force CPU to avoid repeated OOM if we fell through from GPU failure
    )
    gpu_id = gpu_ids[0] if gpu_ids else 0

    t_gpu = time.perf_counter()
    # Process channels sequentially, freeing each input before the next
    # to keep peak memory at ~1 channel + its output maps.

    # DAPI — always present
    dapi_focus_map, dapi_mean_map = compute_ccfs_map(
        dapi, window_size=window_size, use_gpu=use_gpu, gpu_id=gpu_id
    )
    dapi_lap_var_map = compute_laplacian_variance_map(
        dapi,
        window_size=window_size,
        use_gpu=use_gpu,
        gpu_id=gpu_id,
        lap_sigma=lap_sigma,
    )
    del dapi

    # Boundary (channel 1)
    boundary_focus_map: NDArray[np.float32] | None = None
    boundary_mean_map: NDArray[np.float32] | None = None
    if boundary is not None:
        boundary_focus_map, boundary_mean_map = compute_ccfs_map(
            boundary, window_size=window_size, use_gpu=use_gpu, gpu_id=gpu_id
        )
        del boundary

    # IntRNA (channel 2)
    intrna_focus_map: NDArray[np.float32] | None = None
    intrna_mean_map: NDArray[np.float32] | None = None
    if intrna is not None:
        intrna_focus_map, intrna_mean_map = compute_ccfs_map(
            intrna, window_size=window_size, use_gpu=use_gpu, gpu_id=gpu_id
        )
        del intrna

    backend_label = "GPU" if use_gpu else "CPU"
    logging.info(
        f"  [TIMING] {backend_label} compute (single device, all channels): {time.perf_counter() - t_gpu:.1f}s"
    )

    return {
        "dapi_focus_map": dapi_focus_map,
        "dapi_mean_map": dapi_mean_map,
        "dapi_lap_var_map": dapi_lap_var_map,
        "boundary_focus_map": boundary_focus_map,
        "boundary_mean_map": boundary_mean_map,
        "intrna_focus_map": intrna_focus_map,
        "intrna_mean_map": intrna_mean_map,
    }


# ---------------------------------------------------------------------------
# Down-sampling to per-tile DataFrame
# ---------------------------------------------------------------------------


def _build_roi_grid(
    image_shape: tuple[int, int],
    roi_size: int,
    stride: int,
    tissue_mask: NDArray[np.generic] | None,
    downsample_factor: int = 8,
) -> dict[str, Any]:
    """Build the tile grid and compute tissue coverages (one-time setup).

    Returns a dict with keys: ``x1_arr``, ``x2_arr``, ``y1_arr``, ``y2_arr``,
    ``cx``, ``cy``, ``n_rois``, ``height``, ``width``, ``tissue_coverages``,
    ``is_boundary_roi``, ``roi_coords_arr``.
    """
    height, width = image_shape

    # Build ROI grid — vectorized
    x_starts = np.arange(0, width, stride)
    y_starts = np.arange(0, height, stride)
    yy, xx = np.meshgrid(y_starts, x_starts, indexing="ij")
    x1_all = xx.ravel()
    y1_all = yy.ravel()
    x2_all = np.minimum(x1_all + roi_size, width)
    y2_all = np.minimum(y1_all + roi_size, height)

    valid = (x2_all - x1_all >= roi_size // 2) & (y2_all - y1_all >= roi_size // 2)
    x1_arr = x1_all[valid]
    x2_arr = x2_all[valid]
    y1_arr = y1_all[valid]
    y2_arr = y2_all[valid]

    n_rois = len(x1_arr)
    if n_rois == 0:
        raise ValueError(
            f"No tiles generated. Image: {height}x{width}, "
            f"roi_size={roi_size}, stride={stride}."
        )

    roi_coords_arr = np.stack([x1_arr, x2_arr, y1_arr, y2_arr], axis=1)
    cx = (x1_arr + x2_arr) // 2
    cy = (y1_arr + y2_arr) // 2

    # Tissue coverage
    if tissue_mask is not None:
        binary_mask = (tissue_mask.astype(np.float64) > 0).astype(np.float64)
        mask_h, mask_w = binary_mask.shape
        x1_ds = np.clip(x1_arr // downsample_factor, 0, mask_w)
        x2_ds = np.clip((x2_arr - 1) // downsample_factor + 1, 0, mask_w)
        y1_ds = np.clip(y1_arr // downsample_factor, 0, mask_h)
        y2_ds = np.clip((y2_arr - 1) // downsample_factor + 1, 0, mask_h)

        integral = np.zeros((mask_h + 1, mask_w + 1), dtype=np.float64)
        integral[1:, 1:] = np.cumsum(np.cumsum(binary_mask, axis=0), axis=1)
        block_sums = (
            integral[y2_ds, x2_ds]
            - integral[y1_ds, x2_ds]
            - integral[y2_ds, x1_ds]
            + integral[y1_ds, x1_ds]
        )
        block_w = x2_ds - x1_ds
        block_h = y2_ds - y1_ds
        total_pixels_arr = (block_w * block_h).astype(np.float64)
        total_pixels_arr[total_pixels_arr == 0] = 1.0
        tissue_coverages = block_sums / total_pixels_arr
    else:
        tissue_coverages = np.ones(n_rois, dtype=np.float64)

    half_win = roi_size // 2
    is_boundary_roi = (
        (cx < half_win)
        | (cy < half_win)
        | (cx >= width - half_win)
        | (cy >= height - half_win)
    )

    return {
        "x1_arr": x1_arr,
        "x2_arr": x2_arr,
        "y1_arr": y1_arr,
        "y2_arr": y2_arr,
        "cx": cx,
        "cy": cy,
        "n_rois": n_rois,
        "height": height,
        "width": width,
        "tissue_coverages": tissue_coverages,
        "is_boundary_roi": is_boundary_roi,
        "roi_coords_arr": roi_coords_arr,
    }


def downsample_maps_to_roi_dataframe(
    focus_maps: dict[str, NDArray[np.float32] | None],
    tissue_mask: NDArray[np.generic] | None,
    roi_size: int = 35,
    stride: int | None = None,
    downsample_factor: int = 8,
    min_tissue_coverage: float = 0.0,
    image_shape: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """Down-sample pixel-level focus maps to per-tile scalars.

    The output DataFrame has **exactly** the same columns as the legacy
    ``calculate_roi_focusscore()`` function, ensuring full backward
    compatibility.

    For each tile the focus score is obtained by sampling the **centre pixel**
    of the corresponding pixel-level map. Because ``uniform_filter`` at pixel
    ``(cy, cx)`` computes statistics over the surrounding ``window_size x
    window_size`` window, the centre pixel of a ``roi_size x roi_size`` tile
    yields the exact block-level statistic (for non-edge tiles). Edge tiles may
    differ slightly because ``uniform_filter`` uses ``mode='reflect'`` padding
    while the legacy code clips to the image boundary.

    Args:
        focus_maps: Dictionary returned by :func:`compute_all_focus_maps`.
            Required key: ``dapi_focus_map``. All other keys may be ``None``.
        tissue_mask: Labelled tissue mask at down-sampled resolution (e.g.
            level 3). Pass ``None`` to skip tissue filtering (all tiles get
            ``tissue_coverage=1.0``).
        roi_size: Side length of each square tile in pixels.
        stride: Grid spacing in pixels. Defaults to *roi_size* (non-
            overlapping grid).
        downsample_factor: Scale factor between full-resolution coordinates
            and *tissue_mask* coordinates (default 8 for level 3).
        min_tissue_coverage: Minimum tissue fraction for a tile to be
            included. ``0.0`` means any tissue overlap is sufficient. This
            parameter is recorded but **not** used for filtering — all tiles
            are returned so the caller can filter as needed.
        image_shape: ``(height, width)`` of the full-resolution image. If
            ``None`` the shape is inferred from ``dapi_focus_map``.

    Returns:
        :class:`~pandas.DataFrame` with columns identical to the legacy
        ``calculate_roi_focusscore()`` output.

    Raises:
        ValueError: If ``dapi_focus_map`` is missing from *focus_maps* or the
            Tile grid is empty.
    """
    # ------------------------------------------------------------------
    # Unpack maps
    # ------------------------------------------------------------------
    dapi_focus_map = focus_maps.get("dapi_focus_map")
    if dapi_focus_map is None:
        raise ValueError("focus_maps must contain 'dapi_focus_map'.")

    dapi_mean_map = focus_maps.get("dapi_mean_map")
    dapi_lap_var_map = focus_maps.get("dapi_lap_var_map")
    boundary_focus_map = focus_maps.get("boundary_focus_map")
    boundary_mean_map = focus_maps.get("boundary_mean_map")
    intrna_focus_map = focus_maps.get("intrna_focus_map")
    intrna_mean_map = focus_maps.get("intrna_mean_map")

    has_boundary = boundary_focus_map is not None
    has_intrna = intrna_focus_map is not None

    # ------------------------------------------------------------------
    # Image shape
    # ------------------------------------------------------------------
    if image_shape is not None:
        height, width = image_shape
    else:
        height, width = dapi_focus_map.shape

    if stride is None:
        stride = roi_size

    # ------------------------------------------------------------------
    # Build ROI grid — vectorized (identical logic to legacy code)
    # ------------------------------------------------------------------
    t_grid = time.perf_counter()

    x_starts = np.arange(0, width, stride)
    y_starts = np.arange(0, height, stride)
    yy, xx = np.meshgrid(y_starts, x_starts, indexing="ij")
    x1_all = xx.ravel()
    y1_all = yy.ravel()
    x2_all = np.minimum(x1_all + roi_size, width)
    y2_all = np.minimum(y1_all + roi_size, height)

    # Filter out small edge ROIs (same threshold as legacy code)
    valid = (x2_all - x1_all >= roi_size // 2) & (y2_all - y1_all >= roi_size // 2)
    x1_arr = x1_all[valid]
    x2_arr = x2_all[valid]
    y1_arr = y1_all[valid]
    y2_arr = y2_all[valid]

    n_rois = len(x1_arr)
    if n_rois == 0:
        raise ValueError(
            f"No tiles generated from grid. "
            f"Image: {height}x{width}, roi_size={roi_size}, stride={stride}."
        )

    # Keep a structured array for backward-compatible DataFrame columns
    # Columns: x1, x2, y1, y2
    roi_coords_arr = np.stack([x1_arr, x2_arr, y1_arr, y2_arr], axis=1)

    logging.info(
        f"  [TIMING] Tile grid generation ({n_rois} tiles): {time.perf_counter() - t_grid:.1f}s"
    )

    # ------------------------------------------------------------------
    # Tissue coverage — vectorized via block sums
    # ------------------------------------------------------------------
    t_tissue = time.perf_counter()

    tissue_coverages_arr: NDArray[np.float64]
    if tissue_mask is not None:
        # Convert tissue mask to binary
        binary_mask = (tissue_mask.astype(np.float64) > 0).astype(np.float64)
        mask_h, mask_w = binary_mask.shape

        # Compute downsampled ROI coordinates (vectorized)
        x1_ds = x1_arr // downsample_factor
        x2_ds = (x2_arr - 1) // downsample_factor + 1
        y1_ds = y1_arr // downsample_factor
        y2_ds = (y2_arr - 1) // downsample_factor + 1

        # Clip to mask boundaries
        x1_ds = np.clip(x1_ds, 0, mask_w)
        x2_ds = np.clip(x2_ds, 0, mask_w)
        y1_ds = np.clip(y1_ds, 0, mask_h)
        y2_ds = np.clip(y2_ds, 0, mask_h)

        # Check if all ROI blocks have uniform downsampled size
        block_w = x2_ds - x1_ds
        block_h = y2_ds - y1_ds
        uniform_w = int(block_w[0]) if len(block_w) > 0 else 0
        uniform_h = int(block_h[0]) if len(block_h) > 0 else 0
        all_uniform = bool(
            np.all(block_w == uniform_w) and np.all(block_h == uniform_h)
        )

        if all_uniform and uniform_w > 0 and uniform_h > 0:
            # Fast path: use a 2D integral image (summed-area table)
            # to compute block sums in O(1) per tile
            integral = np.zeros((mask_h + 1, mask_w + 1), dtype=np.float64)
            integral[1:, 1:] = np.cumsum(np.cumsum(binary_mask, axis=0), axis=1)
            # Block sum = integral[y2,x2] - integral[y1,x2] - integral[y2,x1] + integral[y1,x1]
            block_sums = (
                integral[y2_ds, x2_ds]
                - integral[y1_ds, x2_ds]
                - integral[y2_ds, x1_ds]
                + integral[y1_ds, x1_ds]
            )
            total_pixels = uniform_w * uniform_h
            tissue_coverages_arr = block_sums / total_pixels
        else:
            # Fallback: still use integral image but handle variable block sizes
            integral = np.zeros((mask_h + 1, mask_w + 1), dtype=np.float64)
            integral[1:, 1:] = np.cumsum(np.cumsum(binary_mask, axis=0), axis=1)
            block_sums = (
                integral[y2_ds, x2_ds]
                - integral[y1_ds, x2_ds]
                - integral[y2_ds, x1_ds]
                + integral[y1_ds, x1_ds]
            )
            total_pixels_arr = (block_w * block_h).astype(np.float64)
            # Avoid division by zero
            total_pixels_arr[total_pixels_arr == 0] = 1.0
            tissue_coverages_arr = block_sums / total_pixels_arr
    else:
        tissue_coverages_arr = np.ones(n_rois, dtype=np.float64)

    logging.info(
        f"  [TIMING] Tissue coverage computation: {time.perf_counter() - t_tissue:.1f}s"
    )

    # ------------------------------------------------------------------
    # Sample centre pixel for each ROI — vectorized fancy indexing
    # ------------------------------------------------------------------
    t_sample = time.perf_counter()

    # Compute centre coordinates for all ROIs at once
    cx = (x1_arr + x2_arr) // 2
    cy = (y1_arr + y2_arr) // 2

    # DAPI channels (always present)
    dapi_focus_scores = dapi_focus_map[cy, cx].astype(np.float64)
    dapi_intensities = (
        dapi_mean_map[cy, cx].astype(np.float64)
        if dapi_mean_map is not None
        else np.zeros(n_rois, dtype=np.float64)
    )
    dapi_lap_vars = (
        dapi_lap_var_map[cy, cx].astype(np.float64)
        if dapi_lap_var_map is not None
        else np.zeros(n_rois, dtype=np.float64)
    )

    # Boundary channel
    if has_boundary:
        boundary_focus_scores = boundary_focus_map[cy, cx].astype(np.float64)  # type: ignore[index]
        boundary_intensities = boundary_mean_map[cy, cx].astype(np.float64)  # type: ignore[index]
    else:
        boundary_focus_scores = None
        boundary_intensities = None

    # IntRNA channel
    if has_intrna:
        intrna_focus_scores = intrna_focus_map[cy, cx].astype(np.float64)  # type: ignore[index]
        intrna_intensities = intrna_mean_map[cy, cx].astype(np.float64)  # type: ignore[index]
    else:
        intrna_focus_scores = None
        intrna_intensities = None

    logging.info(
        f"  [TIMING] Centre-pixel sampling ({n_rois} tiles): {time.perf_counter() - t_sample:.1f}s"
    )

    # ------------------------------------------------------------------
    # Normalize focus scores with RobustScaler
    # ------------------------------------------------------------------
    t_scaler = time.perf_counter()
    scaler_dapi = RobustScaler()
    dapi_focus_scores_norm = scaler_dapi.fit_transform(
        dapi_focus_scores.reshape(-1, 1)
    ).flatten()

    if has_boundary:
        scaler_boundary = RobustScaler()
        boundary_focus_scores_norm: NDArray[np.float64] | None = (
            scaler_boundary.fit_transform(
                boundary_focus_scores.reshape(-1, 1)  # type: ignore[union-attr]
            ).flatten()
        )
    else:
        boundary_focus_scores_norm = None

    if has_intrna:
        scaler_intrna = RobustScaler()
        intrna_focus_scores_norm: NDArray[np.float64] | None = (
            scaler_intrna.fit_transform(
                intrna_focus_scores.reshape(-1, 1)  # type: ignore[union-attr]
            ).flatten()
        )
    else:
        intrna_focus_scores_norm = None

    logging.info(
        f"  [TIMING] RobustScaler normalization: {time.perf_counter() - t_scaler:.1f}s"
    )

    # ------------------------------------------------------------------
    # Mark boundary tiles (centre within roi_size//2 of image edge)
    # ------------------------------------------------------------------
    half_win = roi_size // 2
    is_boundary_roi = (
        (cx < half_win)
        | (cy < half_win)
        | (cx >= width - half_win)
        | (cy >= height - half_win)
    )

    # ------------------------------------------------------------------
    # Build output DataFrame (column order matches legacy code)
    # ------------------------------------------------------------------
    df_data: dict[str, Any] = {
        "roi_id": np.arange(n_rois),
        "x1": roi_coords_arr[:, 0],
        "x2": roi_coords_arr[:, 1],
        "y1": roi_coords_arr[:, 2],
        "y2": roi_coords_arr[:, 3],
        # DAPI focus scores (duplicated for backward compatibility)
        "focus_score": dapi_focus_scores,
        "focus_score_norm": dapi_focus_scores_norm,
        "dapi_focus_score": dapi_focus_scores,
        "dapi_focus_score_norm": dapi_focus_scores_norm,
        # Laplacian variance (DAPI)
        "dapi_lap_var": dapi_lap_vars,
        # Intensities
        "dapi_intensity": dapi_intensities,
        "raw_intensity": dapi_intensities.copy(),  # backward compatibility
        # Boundary channel
        "boundary_focus_score": (boundary_focus_scores if has_boundary else np.nan),
        "boundary_focus_score_norm": (
            boundary_focus_scores_norm if has_boundary else np.nan
        ),
        "boundary_intensity": (boundary_intensities if has_boundary else np.nan),
        # IntRNA channel
        "intrna_focus_score": (intrna_focus_scores if has_intrna else np.nan),
        "intrna_focus_score_norm": (intrna_focus_scores_norm if has_intrna else np.nan),
        "intrna_intensity": (intrna_intensities if has_intrna else np.nan),
        # Tissue coverage
        "tissue_coverage": tissue_coverages_arr,
        "overlaps_tissue": tissue_coverages_arr > 0.0,
        # Boundary flag: True if ROI centre is near image edge
        "is_boundary_roi": is_boundary_roi,
    }

    return pd.DataFrame(df_data)


def calculate_roi_focusscore_without_laplace(
    xoa_morphology_files,
    roi_size=35,
    stride=None,
    tissue_filter=True,
    min_tissue_coverage=0.0,
    downsample_factor=8,
):
    """
    Calculate cell-independent tile-based focus scores using a regular grid.

    Creates a regular lattice/grid of tiles across the entire slide and calculates
    focus scores for each tile independently of cell locations. Calculates focus scores
    for all available channels (DAPI, Boundary, IntRNA).

    Parameters:
    -----------
    xoa_morphology_files : list
        List of paths to morphology image files
    roi_size : int, optional
        Size of square tile in pixels (default: 35)
    stride : int, optional
        Grid spacing in pixels. If None, uses non-overlapping grid (stride = roi_size)
    tissue_filter : bool, optional
        Enable tissue region filtering (default: True)
    min_tissue_coverage : float, optional
        Minimum fraction of tile that must be tissue (default: 0.0, i.e., any tissue overlap)
    downsample_factor : int, optional
        Downsampling factor for tissue mask (default: 8, for level 3)

    Returns:
    --------
    pandas DataFrame
        DataFrame with columns:
        - roi_id: Unique identifier
        - x1, x2, y1, y2: Tile boundaries (full resolution)
        - focus_score: Raw DAPI focus score (std² / mean) [backward compatibility]
        - focus_score_norm: Normalized DAPI focus score [backward compatibility]
        - dapi_focus_score, dapi_focus_score_norm: DAPI focus scores
        - boundary_focus_score, boundary_focus_score_norm: Boundary focus scores (if available)
        - intrna_focus_score, intrna_focus_score_norm: IntRNA focus scores (if available)
        - dapi_intensity: Mean DAPI intensity per tile
        - boundary_intensity: Mean Boundary intensity per tile (if available)
        - intrna_intensity: Mean IntRNA intensity per tile (if available)
        - raw_intensity: Mean DAPI intensity per tile [backward compatibility]
        - tissue_coverage: Fraction of tile that is tissue (0.0-1.0)
    """

    # Set stride (non-overlapping if not specified)
    if stride is None:
        stride = roi_size

    # Load channels from either multi-channel stack or split single-channel files.
    dapi_image, boundary_image, intrna_image = _load_morphology_channels(
        xoa_morphology_files, level=0
    )
    n_channels = 1 + int(boundary_image is not None) + int(intrna_image is not None)

    height, width = dapi_image.shape

    # Check which channels are available
    has_boundary = n_channels > 1
    has_intrna = n_channels > 2

    if has_boundary:
        logging.info(
            f"  Found {n_channels} channels: DAPI, Boundary"
            + (", IntRNA" if has_intrna else "")
        )
    else:
        logging.info(f"  Found {n_channels} channel(s): DAPI only")

    # Generate tissue mask if filtering enabled
    tissue_mask = None
    if tissue_filter:
        # Load downsampled DAPI for tissue mask generation
        small0_ds, _, _ = _load_morphology_channels(xoa_morphology_files, level=3)
        # Generate tissue mask (simplified version - just need whole_sample)
        t = np.percentile(small0_ds, 60)
        thresh1 = small0_ds < t
        thresh2 = small0_ds >= t
        objects = measure.label(thresh1)
        noborder = clear_border(objects)
        holes = morphology.remove_small_objects(noborder, min_size=1500)
        small_objects = noborder ^ holes
        test_mask = measure.label(thresh2) + small_objects
        sample = test_mask > 1
        tissue_mask = measure.label(sample)
        del small0_ds

    # Create grid of ROI coordinates
    n_x = (width // stride) + (1 if width % stride > 0 else 0)
    n_y = (height // stride) + (1 if height % stride > 0 else 0)

    roi_coords = []
    for y_idx in range(n_y):
        for x_idx in range(n_x):
            x1 = x_idx * stride
            x2 = min(x1 + roi_size, width)
            y1 = y_idx * stride
            y2 = min(y1 + roi_size, height)

            # Skip if ROI is too small (at edges)
            if (x2 - x1) < roi_size // 2 or (y2 - y1) < roi_size // 2:
                continue

            roi_coords.append((x1, x2, y1, y2))

    # Calculate tissue coverage for ALL tiles (no filtering)
    # This ensures we always have tiles to process, even if tissue detection fails
    tissue_coverages = []
    if tissue_filter and tissue_mask is not None:
        for x1, x2, y1, y2 in roi_coords:
            # Scale coordinates to downsampled space
            x1_ds = x1 // downsample_factor
            x2_ds = (x2 - 1) // downsample_factor + 1
            y1_ds = y1 // downsample_factor
            y2_ds = (y2 - 1) // downsample_factor + 1

            # Clip to mask boundaries
            x1_ds = max(0, x1_ds)
            x2_ds = min(tissue_mask.shape[1], x2_ds)
            y1_ds = max(0, y1_ds)
            y2_ds = min(tissue_mask.shape[0], y2_ds)

            # Calculate tissue coverage
            roi_mask = tissue_mask[y1_ds:y2_ds, x1_ds:x2_ds]
            tissue_pixels = np.sum(roi_mask > 0)
            total_pixels = roi_mask.size
            coverage = tissue_pixels / total_pixels if total_pixels > 0 else 0.0
            tissue_coverages.append(coverage)
    else:
        # If tissue filtering is disabled, assume all tiles have full tissue coverage
        tissue_coverages = [1.0] * len(roi_coords)

    # Calculate focus scores for all ROIs and all available channels
    n_rois = len(roi_coords)

    # Basic safety check (should never trigger now, but defensive programming)
    if n_rois == 0:
        raise ValueError(
            f"ERROR: No tiles generated from grid!\n"
            f"  - Image dimensions: {height}x{width} (full resolution)\n"
            f"  - Tile size: {roi_size}px, stride: {stride}px\n"
            f"  - This should not happen. Please check image dimensions and tile parameters.\n"
        )

    # Initialize arrays for all channels
    dapi_focus_scores = np.empty(n_rois, dtype=np.float64)
    dapi_intensities = np.empty(n_rois, dtype=np.float64)

    boundary_focus_scores = np.empty(n_rois, dtype=np.float64) if has_boundary else None
    boundary_intensities = np.empty(n_rois, dtype=np.float64) if has_boundary else None

    intrna_focus_scores = np.empty(n_rois, dtype=np.float64) if has_intrna else None
    intrna_intensities = np.empty(n_rois, dtype=np.float64) if has_intrna else None

    # Calculate focus scores for each ROI
    for i, (x1, x2, y1, y2) in enumerate(roi_coords):
        # DAPI channel
        roi_dapi = dapi_image[y1:y2, x1:x2]
        mean_dapi = np.mean(roi_dapi)
        std_dapi = np.std(roi_dapi)
        dapi_focus_scores[i] = (
            (std_dapi * std_dapi) / mean_dapi if mean_dapi > 0 else 0.0
        )
        dapi_intensities[i] = mean_dapi

        # Boundary channel (if available)
        if has_boundary:
            roi_boundary = boundary_image[y1:y2, x1:x2]
            mean_boundary = np.mean(roi_boundary)
            std_boundary = np.std(roi_boundary)
            boundary_focus_scores[i] = (
                (std_boundary * std_boundary) / mean_boundary
                if mean_boundary > 0
                else 0.0
            )
            boundary_intensities[i] = mean_boundary

        # IntRNA channel (if available)
        if has_intrna:
            roi_intrna = intrna_image[y1:y2, x1:x2]
            mean_intrna = np.mean(roi_intrna)
            std_intrna = np.std(roi_intrna)
            intrna_focus_scores[i] = (
                (std_intrna * std_intrna) / mean_intrna if mean_intrna > 0 else 0.0
            )
            intrna_intensities[i] = mean_intrna

    # Normalize focus scores for each channel separately
    # Safety check: Ensure we have data before normalizing
    if len(dapi_focus_scores) == 0:
        raise ValueError(
            f"ERROR: Cannot normalize focus scores - empty array detected!\n"
            f"  - Number of tiles: {n_rois}\n"
            f"  - This should have been caught earlier. Please report this issue.\n"
        )

    scaler_dapi = RobustScaler()
    dapi_focus_scores_norm = scaler_dapi.fit_transform(
        dapi_focus_scores.reshape(-1, 1)
    ).flatten()

    if has_boundary:
        if len(boundary_focus_scores) == 0:
            raise ValueError(
                "ERROR: Cannot normalize boundary focus scores - empty array detected!"
            )
        scaler_boundary = RobustScaler()
        boundary_focus_scores_norm = scaler_boundary.fit_transform(
            boundary_focus_scores.reshape(-1, 1)
        ).flatten()
    else:
        boundary_focus_scores_norm = None

    if has_intrna:
        if len(intrna_focus_scores) == 0:
            raise ValueError(
                "ERROR: Cannot normalize IntRNA focus scores - empty array detected!"
            )
        scaler_intrna = RobustScaler()
        intrna_focus_scores_norm = scaler_intrna.fit_transform(
            intrna_focus_scores.reshape(-1, 1)
        ).flatten()
    else:
        intrna_focus_scores_norm = None

    # Create output DataFrame
    # roi_id: Simple integer IDs (0, 1, 2, ...) for easy indexing
    df_data = {
        "roi_id": range(n_rois),
        "x1": [coords[0] for coords in roi_coords],
        "x2": [coords[1] for coords in roi_coords],
        "y1": [coords[2] for coords in roi_coords],
        "y2": [coords[3] for coords in roi_coords],
        # DAPI focus scores (also kept as focus_score/focus_score_norm for backward compatibility)
        "focus_score": dapi_focus_scores,
        "focus_score_norm": dapi_focus_scores_norm,
        "dapi_focus_score": dapi_focus_scores,
        "dapi_focus_score_norm": dapi_focus_scores_norm,
        "dapi_intensity": dapi_intensities,
        "raw_intensity": dapi_intensities,  # Backward compatibility
        "tissue_coverage": tissue_coverages if tissue_filter else [1.0] * n_rois,
        "overlaps_tissue": [
            coverage > 0.0
            for coverage in (tissue_coverages if tissue_filter else [1.0] * n_rois)
        ],
    }

    # Add Boundary channel data if available
    if has_boundary:
        df_data["boundary_focus_score"] = boundary_focus_scores
        df_data["boundary_focus_score_norm"] = boundary_focus_scores_norm
        df_data["boundary_intensity"] = boundary_intensities
    else:
        df_data["boundary_focus_score"] = np.nan
        df_data["boundary_focus_score_norm"] = np.nan
        df_data["boundary_intensity"] = np.nan

    # Add IntRNA channel data if available
    if has_intrna:
        df_data["intrna_focus_score"] = intrna_focus_scores
        df_data["intrna_focus_score_norm"] = intrna_focus_scores_norm
        df_data["intrna_intensity"] = intrna_intensities
    else:
        df_data["intrna_focus_score"] = np.nan
        df_data["intrna_focus_score_norm"] = np.nan
        df_data["intrna_intensity"] = np.nan

    df_grid_roi = pd.DataFrame(df_data)

    # Clean up
    del dapi_image
    if has_boundary:
        del boundary_image
    if has_intrna:
        del intrna_image

    return df_grid_roi


def calculate_roi_focusscore(
    xoa_morphology_files,
    roi_size=35,
    stride=None,
    tissue_filter=True,
    min_tissue_coverage=0.0,
    downsample_factor=8,
    use_gpu=False,
    gpu_ids=None,
    return_pixel_maps=False,
    lap_sigma: float = 1.0,
):
    """
    Calculate cell-independent tile-based focus scores using a regular grid.

    Uses efficient whole-image convolution operations (uniform_filter + laplace)
    to produce per-pixel focus maps, then samples the centre pixel of each tile
    for backward-compatible per-tile scalars. GPU acceleration is supported via
    CuPy when ``use_gpu=True``. Multi-GPU parallelism is available via
    ``gpu_ids``.

    Parameters:
    -----------
    xoa_morphology_files : list
        List of paths to morphology image files
    roi_size : int, optional
        Size of square tile in pixels (default: 35)
    stride : int, optional
        Grid spacing in pixels. If None, uses non-overlapping grid (stride = roi_size)
    tissue_filter : bool, optional
        Enable tissue region filtering (default: True)
    min_tissue_coverage : float, optional
        Minimum fraction of tile that must be tissue (default: 0.0)
    downsample_factor : int, optional
        Downsampling factor for tissue mask (default: 8, for level 3)
    use_gpu : bool, optional
        Use CuPy GPU backend for convolutions (default: False)
    gpu_ids : list[int] or None, optional
        List of CUDA device IDs for multi-GPU parallelism. Channels are
        distributed across GPUs. If None and use_gpu=True, uses device 0.
    return_pixel_maps : bool, optional
        If True, also return the pixel-level focus map dict (default: False)

    Returns:
    --------
    pandas DataFrame (or tuple of DataFrame, dict if return_pixel_maps=True)
        DataFrame with columns:
        - roi_id: Unique identifier
        - x1, x2, y1, y2: Tile boundaries (full resolution)
        - focus_score: Raw DAPI focus score (std^2 / mean) [backward compat]
        - focus_score_norm: Normalized DAPI focus score [backward compat]
        - dapi_focus_score, dapi_focus_score_norm: DAPI focus scores
        - dapi_lap_var: Laplacian variance (DAPI) per tile
        - boundary_focus_score, boundary_focus_score_norm: Boundary (if avail)
        - intrna_focus_score, intrna_focus_score_norm: IntRNA (if available)
        - dapi_intensity, boundary_intensity, intrna_intensity: Mean per tile
        - raw_intensity: Mean DAPI intensity per tile [backward compat]
        - tissue_coverage: Fraction of tile that is tissue (0.0-1.0)
        - overlaps_tissue: Boolean, tile overlaps any tissue (>0 coverage)
    """
    if stride is None:
        stride = roi_size

    # Resolve GPU configuration early so we can skip the eager load on GPU path.
    _use_gpu = use_gpu
    if gpu_ids is not None and len(gpu_ids) > 0:
        _use_gpu = True

    # Generate tissue mask if filtering enabled (uses small level-3 image only)
    tissue_mask = None
    if tissue_filter:
        small0_ds, _, _ = _load_morphology_channels(xoa_morphology_files, level=3)
        t = np.percentile(small0_ds, 60)
        thresh1 = small0_ds < t
        thresh2 = small0_ds >= t
        objects = measure.label(thresh1)
        noborder = clear_border(objects)
        holes = morphology.remove_small_objects(noborder, min_size=1500)
        small_objects = noborder ^ holes
        test_mask = measure.label(thresh2) + small_objects
        sample = test_mask > 1
        tissue_mask = measure.label(sample)
        del small0_ds

    # ------------------------------------------------------------------
    # GPU path: tiled multi-GPU processing with lazy image loading.
    # Detects channels from lazy TIFF wrappers to avoid loading the
    # full-resolution image into RAM (~50+ GB on large samples).
    # ------------------------------------------------------------------
    if _use_gpu:
        # Open channels as lazy TIFF page wrappers (no pixel data loaded)
        zarr_channels, img_shape = _open_morphology_lazy(xoa_morphology_files, level=0)
        has_boundary = zarr_channels[1] is not None
        has_intrna = zarr_channels[2] is not None
        n_channels = 1 + int(has_boundary) + int(has_intrna)

        gpu_label = f"gpu_ids={gpu_ids}" if gpu_ids else f"gpu={_use_gpu}"
        logging.info(
            f"  Found {n_channels} channel(s), "
            f"computing per-pixel focus maps (window={roi_size}, {gpu_label}, tiled)..."
        )
        logging.info(
            f"  Image shape: {img_shape[0]}x{img_shape[1]}, {n_channels} channel(s)"
        )

        t_gpu = time.perf_counter()

        # Query GPU memory once for adaptive tile sizing
        gpu_mem = cp.cuda.Device(gpu_ids[0]).mem_info[1]  # total VRAM
        logging.info(f"  GPU VRAM: {gpu_mem / 1024**3:.1f} GB per device")

        # Process each channel with tiled multi-GPU
        dapi_result = _compute_channel_maps_tiled(
            zarr_channels[0],
            img_shape,
            roi_size,
            gpu_ids,
            include_laplacian=True,
            lap_sigma=lap_sigma,
            gpu_mem_bytes=gpu_mem,
        )

        boundary_result = None
        if has_boundary and zarr_channels[1] is not None:
            boundary_result = _compute_channel_maps_tiled(
                zarr_channels[1],
                img_shape,
                roi_size,
                gpu_ids,
                include_laplacian=False,
                lap_sigma=lap_sigma,
                gpu_mem_bytes=gpu_mem,
            )

        intrna_result = None
        if has_intrna and zarr_channels[2] is not None:
            intrna_result = _compute_channel_maps_tiled(
                zarr_channels[2],
                img_shape,
                roi_size,
                gpu_ids,
                include_laplacian=False,
                lap_sigma=lap_sigma,
                gpu_mem_bytes=gpu_mem,
            )

        elapsed = time.perf_counter() - t_gpu
        logging.info(
            f"  [TIMING] GPU compute (tiled multi-GPU, {n_channels} ch): {elapsed:.1f}s"
        )

        focus_maps = {
            "dapi_focus_map": dapi_result["focus_map"],
            "dapi_mean_map": dapi_result["mean_map"],
            "dapi_lap_var_map": dapi_result.get("lap_var_map"),
            "boundary_focus_map": (
                boundary_result["focus_map"] if boundary_result else None
            ),
            "boundary_mean_map": (
                boundary_result["mean_map"] if boundary_result else None
            ),
            "intrna_focus_map": (intrna_result["focus_map"] if intrna_result else None),
            "intrna_mean_map": (intrna_result["mean_map"] if intrna_result else None),
        }
        logging.info("  Per-pixel focus maps computed (tiled).")

        df_grid_roi = downsample_maps_to_roi_dataframe(
            focus_maps,
            tissue_mask=tissue_mask,
            roi_size=roi_size,
            stride=stride,
            downsample_factor=downsample_factor,
            min_tissue_coverage=min_tissue_coverage,
        )

        if return_pixel_maps:
            return df_grid_roi, focus_maps

        # Free full-resolution maps now that per-ROI values have been sampled.
        del focus_maps, dapi_result, boundary_result, intrna_result
        del zarr_channels
        gc.collect()
        return df_grid_roi

    # ------------------------------------------------------------------
    # CPU path: incremental per-channel compute→downsample→free to limit
    # peak memory.  At most ~3 maps + convolution intermediates at once.
    # ------------------------------------------------------------------
    # Eagerly load full-resolution channels (CPU path needs numpy arrays)
    dapi_image, boundary_image, intrna_image = _load_morphology_channels(
        xoa_morphology_files, level=0
    )
    has_boundary = boundary_image is not None
    has_intrna = intrna_image is not None
    n_channels = 1 + int(has_boundary) + int(has_intrna)
    logging.info(
        f"  Found {n_channels} channel(s), "
        f"computing per-pixel focus maps (window={roi_size}, gpu=False)..."
    )

    # Build ROI grid once (cheap — just coordinate arrays)
    image_shape = dapi_image.shape[:2]
    grid = _build_roi_grid(
        image_shape, roi_size, stride, tissue_mask, downsample_factor
    )
    cx, cy = grid["cx"], grid["cy"]
    n_rois = grid["n_rois"]
    logging.info(f"  Tile grid: {n_rois:,} tiles")

    # -- DAPI (always present) --
    dapi_focus_map, dapi_mean_map = compute_ccfs_map(
        dapi_image, window_size=roi_size, use_gpu=False, gpu_id=0
    )
    dapi_lap_var_map = compute_laplacian_variance_map(
        dapi_image,
        window_size=roi_size,
        use_gpu=False,
        gpu_id=0,
        lap_sigma=lap_sigma,
    )
    del dapi_image
    # Sample per-tile scalars
    dapi_focus_scores = dapi_focus_map[cy, cx].astype(np.float64)
    dapi_intensities = dapi_mean_map[cy, cx].astype(np.float64)
    dapi_lap_vars = dapi_lap_var_map[cy, cx].astype(np.float64)
    del dapi_lap_var_map  # not needed downstream
    logging.info("  DAPI maps computed and sampled.")

    # -- Boundary --
    boundary_focus_scores = None
    boundary_intensities = None
    boundary_mean_map = None
    if has_boundary:
        b_focus, boundary_mean_map = compute_ccfs_map(
            boundary_image, window_size=roi_size, use_gpu=False, gpu_id=0
        )
        del boundary_image
        boundary_focus_scores = b_focus[cy, cx].astype(np.float64)
        boundary_intensities = boundary_mean_map[cy, cx].astype(np.float64)
        del b_focus  # boundary_mean_map kept for CCFS
        logging.info("  Boundary maps computed and sampled.")
    else:
        del boundary_image

    # -- IntRNA --
    intrna_focus_scores = None
    intrna_intensities = None
    intrna_mean_map = None
    if has_intrna:
        i_focus, intrna_mean_map = compute_ccfs_map(
            intrna_image, window_size=roi_size, use_gpu=False, gpu_id=0
        )
        del intrna_image
        intrna_focus_scores = i_focus[cy, cx].astype(np.float64)
        intrna_intensities = intrna_mean_map[cy, cx].astype(np.float64)
        del i_focus  # intrna_mean_map kept for CCFS
        logging.info("  IntRNA maps computed and sampled.")
    else:
        del intrna_image

    logging.info("  Per-pixel focus maps computed (incremental CPU path).")

    # -- Normalize focus scores with RobustScaler --
    scaler_dapi = RobustScaler()
    dapi_focus_scores_norm = scaler_dapi.fit_transform(
        dapi_focus_scores.reshape(-1, 1)
    ).flatten()

    boundary_focus_scores_norm = None
    if boundary_focus_scores is not None:
        scaler_b = RobustScaler()
        boundary_focus_scores_norm = scaler_b.fit_transform(
            boundary_focus_scores.reshape(-1, 1)
        ).flatten()

    intrna_focus_scores_norm = None
    if intrna_focus_scores is not None:
        scaler_i = RobustScaler()
        intrna_focus_scores_norm = scaler_i.fit_transform(
            intrna_focus_scores.reshape(-1, 1)
        ).flatten()

    # -- Assemble DataFrame --
    rc = grid["roi_coords_arr"]
    df_data: dict[str, Any] = {
        "roi_id": np.arange(n_rois),
        "x1": rc[:, 0],
        "x2": rc[:, 1],
        "y1": rc[:, 2],
        "y2": rc[:, 3],
        "focus_score": dapi_focus_scores,
        "focus_score_norm": dapi_focus_scores_norm,
        "dapi_focus_score": dapi_focus_scores,
        "dapi_focus_score_norm": dapi_focus_scores_norm,
        "dapi_lap_var": dapi_lap_vars,
        "dapi_intensity": dapi_intensities,
        "raw_intensity": dapi_intensities.copy(),
        "boundary_focus_score": boundary_focus_scores if has_boundary else np.nan,
        "boundary_focus_score_norm": boundary_focus_scores_norm
        if has_boundary
        else np.nan,
        "boundary_intensity": boundary_intensities if has_boundary else np.nan,
        "intrna_focus_score": intrna_focus_scores if has_intrna else np.nan,
        "intrna_focus_score_norm": intrna_focus_scores_norm if has_intrna else np.nan,
        "intrna_intensity": intrna_intensities if has_intrna else np.nan,
        "tissue_coverage": grid["tissue_coverages"],
        "overlaps_tissue": grid["tissue_coverages"] > 0.0,
        "is_boundary_roi": grid["is_boundary_roi"],
    }
    df_grid_roi = pd.DataFrame(df_data)

    if return_pixel_maps:
        focus_maps = {
            "dapi_focus_map": dapi_focus_map,
            "dapi_mean_map": dapi_mean_map,
            "dapi_lap_var_map": None,  # already freed
            "boundary_focus_map": None,  # already freed
            "boundary_mean_map": boundary_mean_map,
            "intrna_focus_map": None,  # already freed
            "intrna_mean_map": intrna_mean_map,
        }
        return df_grid_roi, focus_maps
    return df_grid_roi


def calculate_roi_blur_threshold(
    df_grid_roi,
    intensity_threshold=ROI_INTENSITY_THRESHOLD,
    focus_percentile=ROI_FOCUS_SCORE_PERCENTILE,
):
    """
    Calculate tile blur detection threshold from raw focus scores.

    Option B: Exclude tiles with intensity < intensity_threshold, then calculate
    percentile from remaining tissue tiles. This focuses the threshold on tissue regions.

    The threshold is applied to RAW focus scores (not normalized). Classification
    is: blurred if (focus_score <= threshold) OR (intensity < intensity_threshold)

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        DataFrame with 'dapi_focus_score' (raw scores) and 'dapi_intensity'
    intensity_threshold : float, optional
        Minimum intensity to include tiles in percentile calculation (default: ROI_INTENSITY_THRESHOLD)
    focus_percentile : float, optional
        Percentile of raw scores to use as threshold (default: ROI_FOCUS_SCORE_PERCENTILE)

    Returns:
    --------
    float
        Threshold value in raw score units
    """
    # Get intensity column (handle both naming conventions)
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df_grid_roi.columns else "raw_intensity"
    )
    if intensity_col not in df_grid_roi.columns:
        raise ValueError(
            "Intensity column not found. Expected 'dapi_intensity' or 'raw_intensity'"
        )

    # Get focus score column (handle both naming conventions)
    focus_col = (
        "dapi_focus_score"
        if "dapi_focus_score" in df_grid_roi.columns
        else "focus_score"
    )
    if focus_col not in df_grid_roi.columns:
        raise ValueError(
            "Focus score column not found. Expected 'dapi_focus_score' or 'focus_score'"
        )

    # Exclude low-intensity tiles (background) before calculating percentile
    tissue_rois = df_grid_roi[df_grid_roi[intensity_col] >= intensity_threshold]

    if len(tissue_rois) == 0:
        # Fallback: if no tissue tiles, use all tiles
        logging.warning(
            f"  Warning: No tiles with intensity >= {intensity_threshold}, using all tiles for threshold calculation"
        )
        tissue_rois = df_grid_roi

    raw_scores = tissue_rois[focus_col].values

    # Calculate threshold on raw scores from tissue ROIs only
    threshold_raw = np.percentile(raw_scores, focus_percentile)

    return threshold_raw


def fit_focus_gmm(
    df_grid_roi,
    intensity_threshold: float = ROI_INTENSITY_THRESHOLD,
    focus_col_name: str = "dapi_focus_score",
    n_components: int = 2,
    random_state: int = 0,
):
    """
    Fit a Gaussian Mixture Model (GMM) to tile focus scores from tissue tiles
    (intensity >= intensity_threshold) in raw focus-score space.

    The model learns 2 components that approximately correspond to:
    - Lower-focus (blurred) tiles
    - Higher-focus (in-focus) tiles

    This function does NOT classify tiles directly; it only returns the fitted GMM
    and identifies which component is the "blur" component (lower mean).

    Parameters
    ----------
    df_grid_roi : pandas.DataFrame
        DataFrame with at least:
        - focus_col_name (e.g. 'dapi_focus_score' or 'focus_score')
        - 'dapi_intensity' or 'raw_intensity'
    intensity_threshold : float, optional
        Minimum intensity to consider a tile as tissue (default: ROI_INTENSITY_THRESHOLD).
        Tiles below this are excluded from GMM training.
    focus_col_name : str, optional
        Name of the focus-score column to use (default: 'dapi_focus_score').
        If not present, 'focus_score' will be used.
    n_components : int, optional
        Number of Gaussian components for the GMM (default: 2).
    random_state : int, optional
        Random seed for reproducibility (default: 0).

    Returns
    -------
    gmm : sklearn.mixture.GaussianMixture
        Fitted GMM model on (log1p(focus_score)) of tissue tiles.
    blur_component_idx : int
        Index of the GMM component corresponding to blurred tiles (lower mean).
    """
    # Determine intensity column
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df_grid_roi.columns else "raw_intensity"
    )
    if intensity_col not in df_grid_roi.columns:
        raise ValueError(
            "Intensity column not found. Expected 'dapi_intensity' or 'raw_intensity'"
        )

    # Determine focus column
    if focus_col_name not in df_grid_roi.columns:
        # Fallback to generic 'focus_score'
        if "focus_score" not in df_grid_roi.columns:
            raise ValueError(
                f"Focus score column '{focus_col_name}' not found and 'focus_score' not present either"
            )
        focus_col_name = "focus_score"

    # Select tissue tiles for training
    tissue_rois = df_grid_roi[df_grid_roi[intensity_col] >= intensity_threshold].copy()
    if len(tissue_rois) == 0:
        raise ValueError(
            f"No tissue tiles found with intensity >= {intensity_threshold}. "
            f"Cannot fit GMM. Check intensity thresholds or image quality."
        )

    raw_scores = tissue_rois[focus_col_name].values.astype(np.float64)

    # Log-transform to reduce skewness (handle zeros safely)
    x_log = np.log1p(raw_scores).reshape(-1, 1)

    # Fit GMM
    gmm = GaussianMixture(
        n_components=n_components, covariance_type="full", random_state=random_state
    )
    gmm.fit(x_log)

    # Identify which component corresponds to blurred ROIs:
    # assume lower mean log-focus = blur
    means = gmm.means_.flatten()
    blur_component_idx = int(np.argmin(means))

    logging.info("  GMM focus model fitted on tissue tiles")
    logging.info(f"    Number of tissue tiles used: {len(tissue_rois)}")
    logging.info(f"    Component means (log1p focus): {means}")
    logging.info(f"    Blur component index: {blur_component_idx} (lower mean)")

    return gmm, blur_component_idx


def classify_roi_blur(
    df_grid_roi,
    gmm,
    blur_component_idx: int,
    blur_prob_threshold: float = 0.5,
    intensity_threshold: float = ROI_INTENSITY_THRESHOLD,
    focus_col_name: str = "dapi_focus_score",
):
    """
    Classify each tile as blurred or in-focus using a fitted GMM and an intensity safeguard.

    Rules:
    ------
    - Tiles with intensity < intensity_threshold are always marked as blurred
      (low-intensity / background or globally problematic tissue).
    - For tissue tiles (intensity >= threshold), use the GMM posterior probability
      of belonging to the "blur" component:
        - blur_prob = P(component == blur_component_idx | focus_score)
        - Tile is blurred if blur_prob > blur_prob_threshold.

    The classification is added to df_grid_roi in new columns:
    - 'blur_prob_gmm' : posterior probability of being blurred (NaN for low-intensity tiles)
    - 'is_blurred_gmm' : boolean, final classification combining intensity + GMM
    - 'is_low_intensity' : boolean, intensity < intensity_threshold

    Parameters
    ----------
    df_grid_roi : pandas.DataFrame
        DataFrame with at least:
        - focus_col_name (e.g. 'dapi_focus_score' or 'focus_score')
        - 'dapi_intensity' or 'raw_intensity'
    gmm : sklearn.mixture.GaussianMixture
        Fitted GMM model from fit_focus_gmm()
    blur_component_idx : int
        Index of the GMM component corresponding to blurred tiles.
    blur_prob_threshold : float, optional
        Threshold on posterior blur probability to classify a tile as blurred
        (default: 0.5).
    intensity_threshold : float, optional
        Intensity safeguard: tiles below this are auto-blurred (default: ROI_INTENSITY_THRESHOLD).
    focus_col_name : str, optional
        Name of focus-score column used (default: 'dapi_focus_score').

    Returns
    -------
    df_grid_roi : pandas.DataFrame
        Input DataFrame with new columns:
        - 'blur_prob_gmm'
        - 'is_blurred_gmm'
        - 'is_low_intensity'
    """
    df = df_grid_roi.copy()

    # Intensity column
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df.columns else "raw_intensity"
    )
    if intensity_col not in df.columns:
        raise ValueError(
            "Intensity column not found. Expected 'dapi_intensity' or 'raw_intensity'"
        )

    # Focus column
    if focus_col_name not in df.columns:
        if "focus_score" not in df.columns:
            raise ValueError(
                f"Focus score column '{focus_col_name}' not found and 'focus_score' not present either"
            )
        focus_col_name = "focus_score"

    # Intensity-based low-intensity flag
    is_low_intensity = df[intensity_col] < intensity_threshold
    df["is_low_intensity"] = is_low_intensity

    # Initialize columns
    df["blur_prob_gmm"] = np.nan
    df["is_blurred_gmm"] = False

    # Tissue ROIs: intensity >= threshold
    tissue_mask = ~is_low_intensity
    if tissue_mask.any():
        raw_scores = df.loc[tissue_mask, focus_col_name].values.astype(np.float64)
        x_log = np.log1p(raw_scores).reshape(-1, 1)

        # Posterior probabilities of components
        probs = gmm.predict_proba(x_log)
        blur_prob = probs[:, blur_component_idx]

        df.loc[tissue_mask, "blur_prob_gmm"] = blur_prob
        df.loc[tissue_mask, "is_blurred_gmm"] = blur_prob > blur_prob_threshold

    # Low-intensity tiles are always blurred according to the safeguard
    df.loc[is_low_intensity, "is_blurred_gmm"] = True

    # Summary
    total_rois = len(df)
    n_low_int = int(is_low_intensity.sum())
    n_blur = int(df["is_blurred_gmm"].sum())
    pct_blur = (n_blur / total_rois) * 100 if total_rois > 0 else 0.0

    logging.info("  GMM-based blur classification completed")
    logging.info(f"    Total tiles: {total_rois}")
    logging.info(
        f"    Low-intensity tiles (auto-blurred): {n_low_int} ({n_low_int / total_rois * 100:.1f}%)"
    )
    logging.info(f"    Blurred tiles (GMM + intensity): {n_blur} ({pct_blur:.1f}%)")

    return df


def fit_focus_gmm_2d(
    df_grid_roi,
    intensity_threshold: float = ROI_INTENSITY_THRESHOLD,
    focus_col_name: str = "dapi_focus_score",
    n_components: int = 2,
    random_state: int = 0,
):
    """
    Fit a 2D Gaussian Mixture Model (GMM) to tile focus scores using both
    focus score (std²/mean) and Laplacian variance as features.

    The model uses log1p-transformed features:
    - Feature 1: log1p(dapi_focus_score)
    - Feature 2: log1p(dapi_lap_var)

    This allows the model to use both global contrast (focus_score) and
    high-frequency content (Laplacian variance) to distinguish blurred vs in-focus tiles.

    Parameters
    ----------
    df_grid_roi : pandas.DataFrame
        DataFrame with at least:
        - focus_col_name (e.g. 'dapi_focus_score' or 'focus_score')
        - 'dapi_lap_var' (Laplacian variance)
        - 'dapi_intensity' or 'raw_intensity'
    intensity_threshold : float, optional
        Minimum intensity to consider a tile as tissue (default: ROI_INTENSITY_THRESHOLD).
        Tiles below this are excluded from GMM training.
    focus_col_name : str, optional
        Name of the focus-score column to use (default: 'dapi_focus_score').
        If not present, 'focus_score' will be used.
    n_components : int, optional
        Number of Gaussian components for the GMM (default: 2).
    random_state : int, optional
        Random seed for reproducibility (default: 0).

    Returns
    -------
    gmm : sklearn.mixture.GaussianMixture
        Fitted 2D GMM model on [log1p(focus_score), log1p(lap_var)] of tissue tiles.
    blur_component_idx : int
        Index of the GMM component corresponding to blurred tiles (lower mean on focus_score dimension).
    """
    # Determine intensity column
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df_grid_roi.columns else "raw_intensity"
    )
    if intensity_col not in df_grid_roi.columns:
        raise ValueError(
            "Intensity column not found. Expected 'dapi_intensity' or 'raw_intensity'"
        )

    # Determine focus column
    if focus_col_name not in df_grid_roi.columns:
        # Fallback to generic 'focus_score'
        if "focus_score" not in df_grid_roi.columns:
            raise ValueError(
                f"Focus score column '{focus_col_name}' not found and 'focus_score' not present either"
            )
        focus_col_name = "focus_score"

    # Check for Laplacian variance column
    if "dapi_lap_var" not in df_grid_roi.columns:
        raise ValueError(
            "dapi_lap_var column not found. 2D GMM requires Laplacian variance. "
            "Make sure calculate_roi_focusscore() was used (not calculate_roi_focusscore_without_laplace)."
        )

    # Select tissue tiles for training
    tissue_rois = df_grid_roi[df_grid_roi[intensity_col] >= intensity_threshold].copy()
    if len(tissue_rois) == 0:
        raise ValueError(
            f"No tissue tiles found with intensity >= {intensity_threshold}. "
            f"Cannot fit GMM. Check intensity thresholds or image quality."
        )

    # Get both features
    focus_scores = tissue_rois[focus_col_name].values.astype(np.float64)
    lap_vars = tissue_rois["dapi_lap_var"].values.astype(np.float64)

    # Filter out NaN values
    valid_mask = ~(np.isnan(focus_scores) | np.isnan(lap_vars))
    if valid_mask.sum() == 0:
        raise ValueError("No valid tile data (all NaN) for 2D GMM training")

    focus_scores_valid = focus_scores[valid_mask]
    lap_vars_valid = lap_vars[valid_mask]

    # Log-transform both features (clamp negatives to 0 before log1p)
    x_focus_log = np.log1p(np.maximum(focus_scores_valid, 0.0))
    x_lap_log = np.log1p(np.maximum(lap_vars_valid, 0.0))

    # Combine into 2D feature matrix
    x_2d = np.column_stack([x_focus_log, x_lap_log])

    # Fit 2D GMM
    gmm = GaussianMixture(
        n_components=n_components, covariance_type="full", random_state=random_state
    )
    gmm.fit(x_2d)

    # Identify which component corresponds to blurred ROIs:
    # Use the focus_score dimension (first column) - lower mean = blur
    means_focus = gmm.means_[:, 0]  # First dimension (focus_score)
    blur_component_idx = int(np.argmin(means_focus))

    logging.info("  2D GMM focus model fitted on tissue tiles")
    logging.info(f"    Number of tissue tiles used: {valid_mask.sum()}")
    logging.info("    Component means (log1p focus_score, log1p lap_var):")
    for i, mean in enumerate(gmm.means_):
        logging.info(f"      Component {i}: [{mean[0]:.4f}, {mean[1]:.4f}]")
    logging.info(
        f"    Blur component index: {blur_component_idx} (lower mean on focus_score dimension)"
    )

    return gmm, blur_component_idx


def classify_roi_blur_2d(
    df_grid_roi,
    gmm,
    blur_component_idx: int,
    blur_prob_threshold: float = 0.5,
    intensity_threshold: float = ROI_INTENSITY_THRESHOLD,
    focus_col_name: str = "dapi_focus_score",
):
    """
    Classify each tile as blurred or in-focus using a fitted 2D GMM and an intensity safeguard.

    Uses both focus_score and Laplacian variance as features.

    Rules:
    ------
    - Tiles with intensity < intensity_threshold are always marked as blurred
      (low-intensity / background or globally problematic tissue).
    - Tiles missing dapi_lap_var are marked as blurred (cannot use 2D GMM).
    - For tissue tiles (intensity >= threshold) with valid lap_var, use the 2D GMM posterior probability
      of belonging to the "blur" component:
        - blur_prob = P(component == blur_component_idx | focus_score, lap_var)
        - Tile is blurred if blur_prob > blur_prob_threshold.

    The classification is added to df_grid_roi in new columns:
    - 'blur_prob_gmm_2d' : posterior probability of being blurred (NaN for low-intensity or missing lap_var tiles)
    - 'is_blurred_gmm_2d' : boolean, final classification combining intensity + 2D GMM
    - 'is_low_intensity' : boolean, intensity < intensity_threshold (reused if exists)

    Parameters
    ----------
    df_grid_roi : pandas.DataFrame
        DataFrame with at least:
        - focus_col_name (e.g. 'dapi_focus_score' or 'focus_score')
        - 'dapi_lap_var' (Laplacian variance)
        - 'dapi_intensity' or 'raw_intensity'
    gmm : sklearn.mixture.GaussianMixture
        Fitted 2D GMM model from fit_focus_gmm_2d()
    blur_component_idx : int
        Index of the GMM component corresponding to blurred tiles.
    blur_prob_threshold : float, optional
        Threshold on posterior blur probability to classify a tile as blurred
        (default: 0.5).
    intensity_threshold : float, optional
        Intensity safeguard: tiles below this are auto-blurred (default: ROI_INTENSITY_THRESHOLD).
    focus_col_name : str, optional
        Name of focus-score column used (default: 'dapi_focus_score').

    Returns
    -------
    df_grid_roi : pandas.DataFrame
        Input DataFrame with new columns:
        - 'blur_prob_gmm_2d'
        - 'is_blurred_gmm_2d'
        - 'is_low_intensity' (if not already present)
    """
    df = df_grid_roi.copy()

    # Intensity column
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df.columns else "raw_intensity"
    )
    if intensity_col not in df.columns:
        raise ValueError(
            "Intensity column not found. Expected 'dapi_intensity' or 'raw_intensity'"
        )

    # Focus column
    if focus_col_name not in df.columns:
        if "focus_score" not in df.columns:
            raise ValueError(
                f"Focus score column '{focus_col_name}' not found and 'focus_score' not present either"
            )
        focus_col_name = "focus_score"

    # Check for Laplacian variance
    if "dapi_lap_var" not in df.columns:
        raise ValueError(
            "dapi_lap_var column not found. 2D GMM classification requires Laplacian variance."
        )

    # Intensity-based low-intensity flag (reuse if exists, otherwise create)
    if "is_low_intensity" not in df.columns:
        is_low_intensity = df[intensity_col] < intensity_threshold
        df["is_low_intensity"] = is_low_intensity
    else:
        is_low_intensity = df["is_low_intensity"]

    # Initialize columns
    df["blur_prob_gmm_2d"] = np.nan
    df["is_blurred_gmm_2d"] = False

    # Tissue ROIs: intensity >= threshold AND have valid lap_var
    tissue_mask = ~is_low_intensity
    has_lap_var = df["dapi_lap_var"].notna()
    valid_mask = tissue_mask & has_lap_var

    if valid_mask.any():
        # Get both features for valid ROIs
        focus_scores = df.loc[valid_mask, focus_col_name].values.astype(np.float64)
        lap_vars = df.loc[valid_mask, "dapi_lap_var"].values.astype(np.float64)

        # Filter out any remaining NaN (shouldn't happen, but safety check)
        valid_data_mask = ~(np.isnan(focus_scores) | np.isnan(lap_vars))
        if valid_data_mask.sum() > 0:
            focus_scores_valid = focus_scores[valid_data_mask]
            lap_vars_valid = lap_vars[valid_data_mask]

            # Log-transform (clamp negatives to 0, matching fit_focus_gmm_2d)
            x_focus_log = np.log1p(np.maximum(focus_scores_valid, 0.0))
            x_lap_log = np.log1p(np.maximum(lap_vars_valid, 0.0))

            # Combine into 2D feature matrix
            x_2d = np.column_stack([x_focus_log, x_lap_log])

            # Posterior probabilities of components
            probs = gmm.predict_proba(x_2d)
            blur_prob = probs[:, blur_component_idx]

            # Map back to original valid_mask indices
            valid_indices = df.index[valid_mask][valid_data_mask]
            df.loc[valid_indices, "blur_prob_gmm_2d"] = blur_prob
            df.loc[valid_indices, "is_blurred_gmm_2d"] = blur_prob > blur_prob_threshold

    # Low-intensity tiles are always blurred
    df.loc[is_low_intensity, "is_blurred_gmm_2d"] = True

    # Tiles without lap_var are also marked as blurred (cannot use 2D GMM)
    missing_lap_var = ~has_lap_var & tissue_mask
    if missing_lap_var.any():
        df.loc[missing_lap_var, "is_blurred_gmm_2d"] = True

    # Summary
    total_rois = len(df)
    n_low_int = int(is_low_intensity.sum())
    n_missing_lap = int(missing_lap_var.sum()) if missing_lap_var.any() else 0
    n_blur = int(df["is_blurred_gmm_2d"].sum())
    pct_blur = (n_blur / total_rois) * 100 if total_rois > 0 else 0.0

    logging.info("  2D GMM-based blur classification completed")
    logging.info(f"    Total tiles: {total_rois}")
    logging.info(
        f"    Low-intensity tiles (auto-blurred): {n_low_int} ({n_low_int / total_rois * 100:.1f}%)"
    )
    if n_missing_lap > 0:
        logging.info(
            f"    Tiles missing lap_var (auto-blurred): {n_missing_lap} ({n_missing_lap / total_rois * 100:.1f}%)"
        )
    logging.info(f"    Blurred tiles (2D GMM + intensity): {n_blur} ({pct_blur:.1f}%)")

    return df


_MAX_SCATTER_POINTS = 10_000


def _subsample_idx(mask, max_points=None, rng_seed=42):
    """Return indices where *mask* is True, randomly subsampled to *max_points*."""
    if max_points is None:
        max_points = _MAX_SCATTER_POINTS
    idx = np.where(mask)[0]
    if max_points <= 0 or len(idx) <= max_points:
        return idx
    rng = np.random.default_rng(rng_seed)
    return rng.choice(idx, size=max_points, replace=False)


def plot_grid_roi_focus_heatmap(
    df_grid_roi,
    small0,
    figures_dir,
    figures_source_dir,
    threshold=-1.0,
    focus_maps=None,
):
    """
    Create heatmap of grid tile focus scores across the whole tissue.

    When pixel-level ``focus_maps`` are provided, renders smooth per-pixel
    heatmaps via imshow (much faster and higher resolution than Rectangle
    patches). Falls back to the legacy Rectangle-patch approach otherwise.

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: roi_id, x1, x2, y1, y2, focus_score,
        focus_score_norm, raw_intensity, tissue_coverage
    small0 : numpy.ndarray
        Downsampled DAPI image (level 3, 8x downsampling)
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    threshold : float, optional
        Threshold for normalized focus score (default: -1.0)
    focus_maps : dict or None, optional
        Pixel-level focus maps from compute_all_focus_maps(). If provided,
        uses imshow for smooth rendering instead of Rectangle patches.
    """
    from skimage.transform import downscale_local_mean

    downsample_factor = 8
    img_height, img_width = small0.shape
    img_aspect = img_height / img_width if img_width > 0 else 1.0
    panel_width = 7
    panel_height = max(4, panel_width * img_aspect)
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_width + 2, panel_height))

    # Plot 1: Focus score heatmap
    ax = axes[0]
    ax.imshow(
        small0,
        cmap="Greys_r",
        vmax=np.percentile(small0, 99),
        alpha=0.5,
        aspect="auto",
        extent=[0, img_width, img_height, 0],
        origin="upper",
    )

    if focus_maps is not None and focus_maps.get("dapi_focus_map") is not None:
        # Pixel-level heatmap via imshow (fast, smooth)
        dapi_focus = focus_maps["dapi_focus_map"]
        # Downsample pixel map to match small0 resolution
        focus_ds = downscale_local_mean(
            dapi_focus, (downsample_factor, downsample_factor)
        )
        # Clip to small0 dimensions (in case of rounding)
        focus_ds = focus_ds[:img_height, :img_width]
        vmin = np.percentile(focus_ds[focus_ds > 0], 1) if np.any(focus_ds > 0) else 0
        vmax = np.percentile(focus_ds[focus_ds > 0], 99) if np.any(focus_ds > 0) else 1
        ax.imshow(
            focus_ds,
            cmap="viridis",
            alpha=0.6,
            aspect="auto",
            extent=[0, img_width, img_height, 0],
            origin="upper",
            vmin=vmin,
            vmax=vmax,
        )
        sm = plt.cm.ScalarMappable(
            cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax)
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Focus Score (var/mean)", fontsize=12)
        ax.set_title("Per-Pixel Focus Score Heatmap (DAPI)", fontsize=14)
    else:
        # Legacy Rectangle-patch approach
        from matplotlib.patches import Rectangle

        for _, roi in df_grid_roi.iterrows():
            x1_ds = roi["x1"] / downsample_factor
            x2_ds = roi["x2"] / downsample_factor
            y1_ds = roi["y1"] / downsample_factor
            y2_ds = roi["y2"] / downsample_factor
            w = x2_ds - x1_ds
            h = y2_ds - y1_ds
            if w <= 0 or h <= 0:
                continue
            if x1_ds < 0 or y1_ds < 0 or x2_ds > img_width or y2_ds > img_height:
                continue
            norm_score = np.clip((roi["focus_score_norm"] + 3) / 6, 0, 1)
            color = plt.cm.viridis(norm_score)
            rect = Rectangle(
                (x1_ds, y1_ds), w, h, facecolor=color, alpha=0.6, edgecolor="none"
            )
            ax.add_patch(rect)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=-3, vmax=3))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Focus Score (Normalized)", fontsize=12)
        ax.set_title("Grid Tile Focus Score Heatmap (Normalized)", fontsize=14)

    ax.set_xlim(0, img_width)
    ax.set_ylim(img_height, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    # Plot 2: GMM 2D classification (blurred vs in-focus)
    ax = axes[1]
    ax.imshow(
        small0,
        cmap="Greys_r",
        vmax=np.percentile(small0, 99),
        alpha=0.5,
        aspect="auto",
        extent=[0, img_width, img_height, 0],
        origin="upper",
    )

    has_gmm_2d = "is_blurred_gmm_2d" in df_grid_roi.columns

    if (
        focus_maps is not None
        and focus_maps.get("dapi_focus_map") is not None
        and has_gmm_2d
    ):
        # Build blur mask directly at downsampled resolution using vectorized operations
        blur_mask_ds = np.zeros((img_height, img_width), dtype=np.float32)
        blurred_rois = df_grid_roi[df_grid_roi["is_blurred_gmm_2d"]]
        x1_arr = (blurred_rois["x1"].values // downsample_factor).astype(int)
        x2_arr = np.minimum(
            blurred_rois["x2"].values // downsample_factor, img_width
        ).astype(int)
        y1_arr = (blurred_rois["y1"].values // downsample_factor).astype(int)
        y2_arr = np.minimum(
            blurred_rois["y2"].values // downsample_factor, img_height
        ).astype(int)
        for i in range(len(x1_arr)):
            blur_mask_ds[y1_arr[i] : y2_arr[i], x1_arr[i] : x2_arr[i]] = 1.0
        from matplotlib.colors import ListedColormap

        cmap_blur = ListedColormap(["blue", "red"])
        ax.imshow(
            blur_mask_ds,
            cmap=cmap_blur,
            alpha=0.4,
            aspect="auto",
            extent=[0, img_width, img_height, 0],
            origin="upper",
            vmin=0,
            vmax=1,
        )
        ax.set_title(
            "Focus Classification (2D GMM: Blue=In-Focus, Red=Blurred)", fontsize=14
        )
    else:
        # Legacy Rectangle-patch approach
        from matplotlib.patches import Rectangle

        for _, roi in df_grid_roi.iterrows():
            x1_ds = roi["x1"] / downsample_factor
            x2_ds = roi["x2"] / downsample_factor
            y1_ds = roi["y1"] / downsample_factor
            y2_ds = roi["y2"] / downsample_factor
            w = x2_ds - x1_ds
            h = y2_ds - y1_ds
            if w <= 0 or h <= 0:
                continue
            if x1_ds < 0 or y1_ds < 0 or x2_ds > img_width or y2_ds > img_height:
                continue
            if has_gmm_2d:
                color = "blue" if not roi["is_blurred_gmm_2d"] else "red"
            else:
                color = "blue" if roi["focus_score_norm"] > threshold else "red"
            rect = Rectangle(
                (x1_ds, y1_ds), w, h, facecolor=color, alpha=0.6, edgecolor="none"
            )
            ax.add_patch(rect)
        if has_gmm_2d:
            ax.set_title(
                "Grid Tile Focus Score (2D GMM: Blue=In-Focus, Red=Blurred)",
                fontsize=14,
            )
        else:
            ax.set_title(f"Grid Tile Focus Score (Threshold={threshold})", fontsize=14)

    ax.set_xlim(0, img_width)
    ax.set_ylim(img_height, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(
        figures_dir / "grid_roi_focus_heatmap.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV
    df_grid_roi_scaled = df_grid_roi.copy()
    df_grid_roi_scaled["x1_ds"] = df_grid_roi_scaled["x1"] / downsample_factor
    df_grid_roi_scaled["x2_ds"] = df_grid_roi_scaled["x2"] / downsample_factor
    df_grid_roi_scaled["y1_ds"] = df_grid_roi_scaled["y1"] / downsample_factor
    df_grid_roi_scaled["y2_ds"] = df_grid_roi_scaled["y2"] / downsample_factor
    df_grid_roi_scaled.to_csv(
        figures_source_dir / "grid_roi_focus_heatmap.csv", index=False
    )


def _rasterize_heatmap_panel(
    background: np.ndarray,
    df: pd.DataFrame,
    value_col: str,
    cmap,
    vmin: float,
    vmax: float,
    downsample_factor: int = 8,
    alpha: float = 0.6,
    log_norm: bool = False,
) -> np.ndarray:
    """Rasterize a per-tile heatmap onto a morphology background.

    Returns an RGB uint8 array.  Vectorized: computes all tile colours in one
    cmap call, then paints them onto the overlay in a single pass.
    """
    img_h, img_w = background.shape[:2]

    # Normalize background to uint8 RGB
    p99 = float(np.percentile(background, 99)) or 1.0
    bg = np.clip(background / p99 * 255, 0, 255).astype(np.uint8)

    vals = df[value_col].values.astype(np.float64)
    x1s = np.clip(df["x1"].values // downsample_factor, 0, img_w).astype(np.int32)
    x2s = np.clip(df["x2"].values // downsample_factor, 0, img_w).astype(np.int32)
    y1s = np.clip(df["y1"].values // downsample_factor, 0, img_h).astype(np.int32)
    y2s = np.clip(df["y2"].values // downsample_factor, 0, img_h).astype(np.int32)

    if isinstance(cmap, str):
        cmap = plt.cm.get_cmap(cmap)

    # Vectorized normalisation
    finite = np.isfinite(vals)
    t = np.zeros_like(vals)
    if log_norm:
        safe = np.where(finite, np.maximum(vals, vmin), vmin)
        log_range = np.log10(vmax) - np.log10(vmin)
        t = (np.log10(safe) - np.log10(vmin)) / log_range if log_range > 0 else t
    else:
        t = (vals - vmin) / (vmax - vmin) if vmax > vmin else t
    t = np.clip(t, 0, 1)

    # Vectorised colormap: one call for all tiles → (N, 4) RGBA uint8
    rgba_all = (cmap(t)[:, :3] * 255).astype(np.uint8)  # (N, 3) RGB

    # Alpha-blend onto background
    result = np.stack([bg, bg, bg], axis=-1)  # (H, W, 3) uint8
    a_uint8 = int(alpha * 255)
    inv_a = 255 - a_uint8

    for i in range(len(vals)):
        if not finite[i]:
            continue
        r1, r2, c1, c2 = y1s[i], y2s[i], x1s[i], x2s[i]
        if r2 <= r1 or c2 <= c1:
            continue
        # In-place uint8 alpha-blend (no float conversion)
        tile = result[r1:r2, c1:c2]
        tile[:] = (
            (tile.astype(np.uint16) * inv_a + rgba_all[i].astype(np.uint16) * a_uint8)
            >> 8
        ).astype(np.uint8)

    return result


def _make_colorbar(
    cmap,
    vmin: float,
    vmax: float,
    height: int,
    width: int = 30,
    label: str = "",
    log_norm: bool = False,
    threshold: float | None = None,
) -> np.ndarray:
    """Render a vertical colorbar as an RGB uint8 array."""
    bar_h = height - 40  # leave room for labels
    bar = np.zeros((bar_h, width, 3), dtype=np.uint8)
    if isinstance(cmap, str):
        cmap = plt.cm.get_cmap(cmap)
    for row in range(bar_h):
        t = 1.0 - row / max(bar_h - 1, 1)  # top=vmax, bottom=vmin
        rgba = cmap(t)
        bar[row, :] = [int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)]

    # Assemble into a wider strip with labels
    total_w = width + 80
    result = np.ones((height, total_w, 3), dtype=np.uint8) * 255
    y_off = 20
    result[y_off : y_off + bar_h, 5 : 5 + width] = bar

    # Draw labels via PIL
    img = Image.fromarray(result)
    draw = ImageDraw.Draw(img)
    draw.text((5 + width + 4, y_off - 2), f"{vmax:.2g}", fill=(0, 0, 0))
    draw.text((5 + width + 4, y_off + bar_h - 12), f"{vmin:.2g}", fill=(0, 0, 0))
    if label:
        # Vertical-ish label at bottom
        draw.text((2, y_off + bar_h + 4), label[:18], fill=(0, 0, 0))
    if threshold is not None and vmax > vmin:
        frac = (threshold - vmin) / (vmax - vmin)
        y_line = y_off + int((1 - frac) * bar_h)
        draw.line([(0, y_line), (5 + width, y_line)], fill=(0, 0, 0), width=2)

    return np.array(img)


def _assemble_two_panel_pil(
    left: np.ndarray,
    right: np.ndarray,
    left_cbar: np.ndarray,
    right_cbar: np.ndarray,
    left_title: str,
    right_title: str,
    gap: int = 20,
) -> Image.Image:
    """Concatenate two panels + colorbars horizontally, add titles."""
    h = max(left.shape[0], right.shape[0])
    title_h = 30
    total_h = h + title_h

    # Resize colorbars to panel height
    left_cbar_img = Image.fromarray(left_cbar).resize(
        (left_cbar.shape[1], h), Image.NEAREST
    )
    right_cbar_img = Image.fromarray(right_cbar).resize(
        (right_cbar.shape[1], h), Image.NEAREST
    )

    total_w = (
        left.shape[1] + left_cbar.shape[1] + gap + right.shape[1] + right_cbar.shape[1]
    )
    canvas = Image.new("RGB", (total_w, total_h), (255, 255, 255))

    x = 0
    canvas.paste(Image.fromarray(left), (x, title_h))
    x += left.shape[1]
    canvas.paste(left_cbar_img, (x, title_h))
    x += left_cbar.shape[1] + gap
    right_x = x
    canvas.paste(Image.fromarray(right), (x, title_h))
    x += right.shape[1]
    canvas.paste(right_cbar_img, (x, title_h))

    draw = ImageDraw.Draw(canvas)
    draw.text((10, 5), left_title, fill=(0, 0, 0))
    draw.text((right_x + 10, 5), right_title, fill=(0, 0, 0))

    return canvas


def plot_grid_roi_blur_heatmap(
    df_grid_roi,
    small0,
    figures_dir,
    figures_source_dir,
    blur_prob_threshold=0.5,
):
    """Spatial blur probability and classification heatmap (PIL rendering).

    Two-panel figure:
      Left:  GMM blur probability (continuous 0-1), coloured RdYlGn_r so red = blurred.
      Right: Binary blur classification (in-focus vs blurred) from 2D GMM.

    Saved as ``figures/grid_roi_blur_heatmap.png``.
    """
    prob_col = (
        "blur_prob_gmm_2d"
        if "blur_prob_gmm_2d" in df_grid_roi.columns
        else "blur_prob_gmm"
        if "blur_prob_gmm" in df_grid_roi.columns
        else None
    )
    class_col = (
        "is_blurred_gmm_2d"
        if "is_blurred_gmm_2d" in df_grid_roi.columns
        else "is_blurred_gmm"
        if "is_blurred_gmm" in df_grid_roi.columns
        else None
    )
    if prob_col is None and class_col is None:
        logging.info(
            "Skipping blur heatmap: no blur_prob or is_blurred columns in df_grid_roi"
        )
        return

    # ── Left panel: blur probability (continuous) ──
    if prob_col is not None:
        left = _rasterize_heatmap_panel(
            small0, df_grid_roi, prob_col, "RdYlGn_r", 0.0, 1.0
        )
        left_cbar = _make_colorbar(
            "RdYlGn_r",
            0.0,
            1.0,
            left.shape[0],
            label="Blur Prob (GMM)",
            threshold=blur_prob_threshold,
        )
        left_title = "Blur Probability per Tile"
    else:
        h, w = small0.shape
        left = np.ones((h, w, 3), dtype=np.uint8) * 200
        left_cbar = np.ones((h, 110, 3), dtype=np.uint8) * 255
        left_title = "Blur probability not available"

    # ── Right panel: binary blur classification ──
    if class_col is not None:
        from matplotlib.colors import ListedColormap

        binary_cmap = ListedColormap(["#2196F3", "#F44336"])
        right = _rasterize_heatmap_panel(
            small0, df_grid_roi, class_col, binary_cmap, 0.0, 1.0
        )
        right_cbar = _make_colorbar(
            binary_cmap, 0.0, 1.0, right.shape[0], label="Focus / Blur"
        )
        right_title = "Blur Classification (GMM)"
    else:
        h, w = small0.shape
        right = np.ones((h, w, 3), dtype=np.uint8) * 200
        right_cbar = np.ones((h, 110, 3), dtype=np.uint8) * 255
        right_title = "Blur classification not available"

    img = _assemble_two_panel_pil(
        left, right, left_cbar, right_cbar, left_title, right_title
    )
    img.save(figures_dir / "grid_roi_blur_heatmap.png")

    # Save source data
    cols = ["roi_id", "x1", "x2", "y1", "y2"]
    if prob_col:
        cols.append(prob_col)
    if class_col:
        cols.append(class_col)
    if "tissue_coverage" in df_grid_roi.columns:
        cols.append("tissue_coverage")
    export = df_grid_roi[[c for c in cols if c in df_grid_roi.columns]].copy()
    export.to_csv(figures_source_dir / "grid_roi_blur_heatmap.csv", index=False)


def plot_snr_roi_heatmap(
    df_grid_roi,
    small0,
    figures_dir,
    figures_source_dir,
):
    """Paint per-tile transcript SNR spatially on the DAPI background (PIL rendering).

    Two-panel figure:
      Left  – neg_pct (fraction of negative-control transcripts per tile)
      Right – roi_tx_snr_ratio (real / negative transcript ratio, log scale)
    """
    # Filter to tiles with transcripts assigned
    df_plot = df_grid_roi[df_grid_roi["snr_total_tx"] > 0].copy()
    if df_plot.empty:
        logging.warning("No tiles with transcripts — skipping SNR heatmap.")
        return

    downsample_factor = 8

    # ── Left panel: neg_pct ──
    neg_vals = df_plot["neg_pct"].values
    neg_vmax = max(0.30, float(np.nanpercentile(neg_vals, 99)))

    left = _rasterize_heatmap_panel(
        small0, df_plot, "neg_pct", "RdYlGn_r", 0.0, neg_vmax
    )
    left_cbar = _make_colorbar(
        "RdYlGn_r", 0.0, neg_vmax, left.shape[0], label="neg_pct (frac)"
    )

    # ── Right panel: roi_tx_snr_ratio (log scale) ──
    ratio_vals = df_plot["roi_tx_snr_ratio"].values
    finite_ratios = ratio_vals[np.isfinite(ratio_vals) & (ratio_vals > 0)]
    ratio_vmax = (
        max(100.0, float(np.nanpercentile(finite_ratios, 99)))
        if len(finite_ratios)
        else 100.0
    )

    right = _rasterize_heatmap_panel(
        small0, df_plot, "roi_tx_snr_ratio", "RdYlGn", 1.0, ratio_vmax, log_norm=True
    )
    right_cbar = _make_colorbar(
        "RdYlGn",
        1.0,
        ratio_vmax,
        right.shape[0],
        label="SNR ratio",
        log_norm=True,
    )

    img = _assemble_two_panel_pil(
        left,
        right,
        left_cbar,
        right_cbar,
        "Negative Probe Fraction per Tile",
        "Transcript SNR Ratio per Tile (log scale)",
    )
    img.save(figures_dir / "snr_heatmap.png")

    # Save source data
    _src_cols = [
        "roi_id",
        "x1",
        "x2",
        "y1",
        "y2",
        "neg_pct",
        "roi_tx_snr_ratio",
        "roi_tx_snr_log",
        "tissue_coverage",
    ]
    df_src = df_plot[[c for c in _src_cols if c in df_plot.columns]].copy()
    df_src["x1_ds"] = df_src["x1"] / downsample_factor
    df_src["x2_ds"] = df_src["x2"] / downsample_factor
    df_src["y1_ds"] = df_src["y1"] / downsample_factor
    df_src["y2_ds"] = df_src["y2"] / downsample_factor
    df_src.to_csv(figures_source_dir / "snr_heatmap.csv", index=False)


def plot_cross_section_concordance(
    df_grid_roi,
    figures_dir,
    figures_source_dir,
    snr_thresholds=None,
):
    """Cross-section concordance: image quality vs transcript SNR per tile.

    Two-panel scatter showing how different quality lenses relate:
      Left  – focus_score vs roi_tx_snr_ratio (optical quality → transcript quality)
      Right – dapi_intensity vs neg_pct (signal strength → noise contamination)

    Only tissue tiles with transcripts are included.
    """
    required = {"focus_score", "roi_tx_snr_ratio", "neg_pct", "dapi_intensity"}
    missing = required - set(df_grid_roi.columns)
    if missing:
        logging.warning(
            "Skipping cross-section concordance plot: missing columns %s", missing
        )
        return

    # Filter to tissue ROIs with transcripts
    df = df_grid_roi.copy()
    if "overlaps_tissue" in df.columns:
        df = df[df["overlaps_tissue"]]
    if "snr_total_tx" in df.columns:
        df = df[df["snr_total_tx"] > 0]
    mask = (
        df["focus_score"].notna()
        & df["roi_tx_snr_ratio"].notna()
        & df["neg_pct"].notna()
        & df["dapi_intensity"].notna()
    )
    df = df[mask]
    if len(df) < 10:
        logging.warning(
            "Skipping cross-section concordance: too few valid tiles (%d)", len(df)
        )
        return

    from scipy.stats import spearmanr

    _t = snr_thresholds or {}

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # --- Left panel: focus_score vs roi_tx_snr_ratio ---
    ax = axes[0]
    focus_vals = df["focus_score"].values
    snr_vals = df["roi_tx_snr_ratio"].values

    sidx = _subsample_idx(np.ones(len(df), dtype=bool))
    # Color by GMM 2D classification if available
    if "is_blurred_gmm_2d" in df.columns and df["is_blurred_gmm_2d"].notna().any():
        colors = np.where(
            df["is_blurred_gmm_2d"].values[sidx].astype(bool), "#E57373", "#64B5F6"
        )
    else:
        colors = "#64B5F6"

    ax.scatter(
        np.log1p(focus_vals[sidx]),
        np.log1p(snr_vals[sidx]),
        s=8,
        alpha=0.4,
        c=colors,
        edgecolors="none",
        rasterized=True,
    )
    rho_fs, pval_fs = spearmanr(focus_vals, snr_vals)
    ax.set_xlabel("log₁₊(Focus Score)", fontsize=12)
    ax.set_ylabel("log₁₊(Transcript SNR Ratio)", fontsize=12)
    ax.set_title("Optical Quality vs Transcript Quality", fontsize=13)
    ax.text(
        0.05,
        0.95,
        f"Spearman ρ = {rho_fs:.3f}\nn = {len(df):,} tiles",
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
    )
    # Threshold lines
    ratio_warn = float(_t.get("ratio_warn", 3.0))
    ratio_fail = float(_t.get("ratio_fail", 1.5))
    ax.axhline(
        np.log1p(ratio_warn),
        ls="--",
        lw=1,
        color="orange",
        alpha=0.8,
        label=f"WARN = {ratio_warn}",
    )
    ax.axhline(
        np.log1p(ratio_fail),
        ls="--",
        lw=1,
        color="red",
        alpha=0.8,
        label=f"FAIL = {ratio_fail}",
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, ls="--", alpha=0.3)

    # --- Right panel: dapi_intensity vs neg_pct ---
    ax = axes[1]
    int_vals = df["dapi_intensity"].values
    neg_vals = df["neg_pct"].values

    if "is_blurred_gmm_2d" in df.columns and df["is_blurred_gmm_2d"].notna().any():
        colors_r = np.where(
            df["is_blurred_gmm_2d"].values[sidx].astype(bool), "#E57373", "#64B5F6"
        )
    else:
        colors_r = "#64B5F6"

    ax.scatter(
        np.log1p(int_vals[sidx]),
        neg_vals[sidx],
        s=8,
        alpha=0.4,
        c=colors_r,
        edgecolors="none",
        rasterized=True,
    )
    rho_in, pval_in = spearmanr(int_vals, neg_vals)
    ax.set_xlabel("log₁₊(DAPI Intensity)", fontsize=12)
    ax.set_ylabel("Negative Probe Fraction", fontsize=12)
    ax.set_title("Signal Strength vs Noise Contamination", fontsize=13)
    ax.text(
        0.95,
        0.95,
        f"Spearman ρ = {rho_in:.3f}\nn = {len(df):,} tiles",
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
    )
    neg_warn = float(_t.get("neg_pct_warn", 0.15))
    neg_fail = float(_t.get("neg_pct_fail", 0.30))
    ax.axhline(
        neg_warn, ls="--", lw=1, color="orange", alpha=0.8, label=f"WARN = {neg_warn}"
    )
    ax.axhline(
        neg_fail, ls="--", lw=1, color="red", alpha=0.8, label=f"FAIL = {neg_fail}"
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, ls="--", alpha=0.3)

    # Add GMM legend if coloured
    if "is_blurred_gmm_2d" in df.columns and df["is_blurred_gmm_2d"].notna().any():
        from matplotlib.patches import Patch

        for a in axes:
            handles = a.get_legend_handles_labels()[0]
            handles.extend(
                [
                    Patch(facecolor="#64B5F6", label="In-Focus (GMM)"),
                    Patch(facecolor="#E57373", label="Blurred (GMM)"),
                ]
            )
            a.legend(
                handles=handles,
                fontsize=8,
                loc="lower right" if a is axes[0] else "upper right",
            )

    plt.tight_layout()
    plt.savefig(
        figures_dir / "cross_section_concordance.png", dpi=200, bbox_inches="tight"
    )
    plt.savefig(
        figures_dir / "cross_section_concordance.pdf", dpi=200, bbox_inches="tight"
    )
    plt.close(fig)

    # Source data
    src_cols = [
        "roi_id",
        "focus_score",
        "dapi_intensity",
        "roi_tx_snr_ratio",
        "neg_pct",
    ]
    if "is_blurred_gmm_2d" in df.columns:
        src_cols.append("is_blurred_gmm_2d")
    if "tissue_coverage" in df.columns:
        src_cols.append("tissue_coverage")
    df[[c for c in src_cols if c in df.columns]].to_csv(
        figures_source_dir / "cross_section_concordance.csv", index=False
    )
    logging.info(
        "Cross-section concordance: focus-vs-SNR ρ=%.3f, intensity-vs-neg ρ=%.3f (%d tiles)",
        rho_fs,
        rho_in,
        len(df),
    )


def plot_roi_focus_vs_intensity(
    df_grid_roi, figures_dir, figures_source_dir, threshold=-1.0
):
    """
    Create scatter plot showing focus score vs raw DAPI intensity for tile threshold analysis.
    Uses log scale for intensity to better visualize the relationship.

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: focus_score, focus_score_norm, raw_intensity (or dapi_intensity)
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    threshold : float, optional
        Threshold for normalized focus score (default: -1.0)
    """
    # Use dapi_intensity if available, otherwise raw_intensity
    intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df_grid_roi.columns else "raw_intensity"
    )

    # Calculate log10 of intensity (add small epsilon to avoid log(0))
    intensity_values = df_grid_roi[intensity_col].values
    log_intensity = np.log10(
        intensity_values + 1e-10
    )  # Add small epsilon to handle any zeros

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: Raw focus score vs log intensity
    ax = axes[0]
    # Filter out invalid values
    valid_mask = (
        np.isfinite(log_intensity)
        & np.isfinite(df_grid_roi["focus_score"].values)
        & np.isfinite(df_grid_roi["focus_score_norm"].values)
    )
    if valid_mask.sum() > 0:
        sidx = _subsample_idx(valid_mask)
        scatter = ax.scatter(
            log_intensity[sidx],
            df_grid_roi["focus_score"].values[sidx],
            alpha=0.5,
            s=10,
            c=df_grid_roi["focus_score_norm"].values[sidx],
            cmap="viridis",
            edgecolors="none",
            rasterized=True,
        )
        plt.colorbar(scatter, ax=ax, label="Normalized Focus Score")
        # Ensure axes are visible
        ax.set_xlim(
            log_intensity[valid_mask].min() - 0.1, log_intensity[valid_mask].max() + 0.1
        )
        ax.set_ylim(
            df_grid_roi.loc[valid_mask, "focus_score"].min() * 0.9,
            df_grid_roi.loc[valid_mask, "focus_score"].max() * 1.1,
        )
    else:
        logging.warning("  Warning: No valid data points for scatter plot")
        ax.text(
            0.5,
            0.5,
            "No valid data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
    ax.set_xlabel("Log10(Raw DAPI Intensity) (mean per tile)", fontsize=12)
    ax.set_ylabel("Focus Score (Raw: std²/mean)", fontsize=12)
    ax.set_title("Focus Score vs Log10(Raw DAPI Intensity)", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.7)

    # Plot 2: Normalized focus score vs log intensity (with GMM 2D classification)
    ax = axes[1]
    # Filter out invalid values
    valid_mask_plot2 = np.isfinite(log_intensity) & np.isfinite(
        df_grid_roi["focus_score_norm"].values
    )

    if valid_mask_plot2.sum() == 0:
        logging.warning("  Warning: No valid data points for plot 2")
        ax.text(
            0.5,
            0.5,
            "No valid data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
    else:
        # Color by GMM 2D classification if available, otherwise fall back to threshold
        has_gmm_2d = "is_blurred_gmm_2d" in df_grid_roi.columns

        if has_gmm_2d:
            # Color by GMM 2D classification
            valid_gmm_mask = valid_mask_plot2 & df_grid_roi["is_blurred_gmm_2d"].notna()
            if valid_gmm_mask.sum() > 0:
                sidx = _subsample_idx(valid_gmm_mask)
                is_blurred_sub = (
                    df_grid_roi["is_blurred_gmm_2d"].values[sidx].astype(bool)
                )
                colors = np.where(is_blurred_sub, "red", "blue")
                ax.scatter(
                    log_intensity[sidx],
                    df_grid_roi["focus_score_norm"].values[sidx],
                    alpha=0.5,
                    s=10,
                    c=colors,
                    edgecolors="none",
                    rasterized=True,
                )

                # Set axes limits from full data to ensure visibility
                ax.set_xlim(
                    log_intensity[valid_gmm_mask].min() - 0.1,
                    log_intensity[valid_gmm_mask].max() + 0.1,
                )
                ax.set_ylim(
                    df_grid_roi.loc[valid_gmm_mask, "focus_score_norm"].min() - 0.5,
                    df_grid_roi.loc[valid_gmm_mask, "focus_score_norm"].max() + 0.5,
                )

                # Add text annotation for GMM 2D (counts from full data, not subsample)
                n_blurred = df_grid_roi.loc[valid_gmm_mask, "is_blurred_gmm_2d"].sum()
                n_in_focus = (
                    ~df_grid_roi.loc[valid_gmm_mask, "is_blurred_gmm_2d"]
                ).sum()
                pct_blurred = (
                    n_blurred / valid_gmm_mask.sum() * 100
                    if valid_gmm_mask.sum() > 0
                    else 0
                )
                ax.text(
                    0.05,
                    0.95,
                    f"2D GMM Classification\nRed: Blurred ({n_blurred:,}, {pct_blurred:.1f}%)\nBlue: In-Focus ({n_in_focus:,}, {100 - pct_blurred:.1f}%)",
                    transform=ax.transAxes,
                    fontsize=11,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
                )
                ax.set_title(
                    "Normalized Focus Score vs Log10(Raw DAPI Intensity) (2D GMM Classification)",
                    fontsize=14,
                )
        else:
            # Fallback to threshold method
            sidx = _subsample_idx(valid_mask_plot2)
            scores_sub = df_grid_roi["focus_score_norm"].values[sidx]
            colors = np.where(scores_sub <= threshold, "red", "blue")
            ax.scatter(
                log_intensity[sidx],
                scores_sub,
                alpha=0.5,
                s=10,
                c=colors,
                edgecolors="none",
                rasterized=True,
            )

            # Set axes limits from full data to ensure visibility
            ax.set_xlim(
                log_intensity[valid_mask_plot2].min() - 0.1,
                log_intensity[valid_mask_plot2].max() + 0.1,
            )
            ax.set_ylim(
                df_grid_roi.loc[valid_mask_plot2, "focus_score_norm"].min() - 0.5,
                df_grid_roi.loc[valid_mask_plot2, "focus_score_norm"].max() + 0.5,
            )

            # Add threshold line
            ax.axhline(
                y=threshold,
                color="black",
                linestyle="--",
                linewidth=2,
                label=f"Threshold ({threshold})",
            )

            # Add text annotation for threshold
            ax.text(
                0.05,
                0.95,
                f"Threshold: {threshold}\nRed: Blurred (≤{threshold})\nBlue: In-Focus (>{threshold})",
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
            )
            ax.set_title(
                "Normalized Focus Score vs Log10(Raw DAPI Intensity) (Thresholded)",
                fontsize=14,
            )

    ax.set_xlabel("Log10(Raw DAPI Intensity) (mean per tile)", fontsize=12)
    ax.set_ylabel("Focus Score (Normalized)", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(
        figures_dir / "roi_focus_vs_intensity.pdf", dpi=300, bbox_inches="tight"
    )
    plt.savefig(
        figures_dir / "roi_focus_vs_intensity.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV (include both raw and log intensity)
    df_scatter = df_grid_roi[
        ["focus_score", "focus_score_norm", intensity_col, "tissue_coverage"]
    ].copy()
    df_scatter["log10_intensity"] = log_intensity
    df_scatter["is_low_nuclear_texture"] = df_scatter["focus_score_norm"] <= threshold
    df_scatter.to_csv(figures_source_dir / "roi_focus_vs_intensity.csv", index=False)

    # Print correlation statistics (using log intensity)
    correlation = np.corrcoef(log_intensity, df_grid_roi["focus_score_norm"])[0, 1]
    logging.info(
        f"  Correlation (log10(intensity) vs normalized focus score): {correlation:.4f}"
    )


def plot_roi_focus_distribution(
    df_grid_roi, figures_dir, figures_source_dir, threshold=-1.0
):
    """
    Create histogram and density plot of focus score values across lattice squares.

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: focus_score, focus_score_norm
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    threshold : float, optional
        Threshold for normalized focus score (default: -1.0)
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Histogram of raw focus scores
    ax = axes[0, 0]
    ax.hist(
        df_grid_roi["focus_score"],
        bins=50,
        alpha=0.7,
        edgecolor="black",
        color="steelblue",
    )
    ax.set_xlabel("Focus Score (Raw: std²/mean)", fontsize=12)
    ax.set_ylabel("Number of tiles", fontsize=12)
    ax.set_title("Distribution of Raw Focus Scores", fontsize=14)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # Add statistics text
    mean_fs = df_grid_roi["focus_score"].mean()
    median_fs = df_grid_roi["focus_score"].median()
    stats_text = f"Mean: {mean_fs:.4f}\nMedian: {median_fs:.4f}"
    ax.text(
        0.95,
        0.95,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
        fontsize=10,
    )

    # Plot 2: Histogram of normalized focus scores with GMM 2D classification
    ax = axes[0, 1]
    # Use GMM 2D classification if available, otherwise fall back to threshold
    has_gmm_2d = "is_blurred_gmm_2d" in df_grid_roi.columns

    if has_gmm_2d:
        blurred = df_grid_roi[df_grid_roi["is_blurred_gmm_2d"]]
        in_focus = df_grid_roi[~df_grid_roi["is_blurred_gmm_2d"]]
        title_suffix = "2D GMM Classification"
    else:
        blurred = df_grid_roi[df_grid_roi["focus_score_norm"] <= threshold]
        in_focus = df_grid_roi[df_grid_roi["focus_score_norm"] > threshold]
        title_suffix = f"Threshold: {threshold}"

    ax.hist(
        [blurred["focus_score_norm"], in_focus["focus_score_norm"]],
        bins=50,
        alpha=0.7,
        edgecolor="black",
        color=["red", "blue"],
        label=["Blurred", "In-Focus"],
        stacked=False,
    )

    # Add threshold line only if using threshold method
    if not has_gmm_2d:
        ax.axvline(
            x=threshold,
            color="black",
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({threshold})",
        )

    ax.set_xlabel("Focus Score (Normalized)", fontsize=12)
    ax.set_ylabel("Number of tiles", fontsize=12)
    ax.set_title(
        f"Distribution of Normalized Focus Scores ({title_suffix})", fontsize=14
    )
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # Add statistics text
    mean_norm = df_grid_roi["focus_score_norm"].mean()
    median_norm = df_grid_roi["focus_score_norm"].median()
    pct_blurred = len(blurred) / len(df_grid_roi) * 100
    stats_text = (
        f"Mean: {mean_norm:.4f}\nMedian: {median_norm:.4f}\nBlurred: {pct_blurred:.1f}%"
    )
    ax.text(
        0.95,
        0.95,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
        fontsize=10,
    )

    # Plot 3: Density plot of raw focus scores
    ax = axes[1, 0]
    ax.hist(
        df_grid_roi["focus_score"],
        bins=50,
        alpha=0.7,
        edgecolor="black",
        color="steelblue",
        density=True,
    )
    ax.set_xlabel("Focus Score (Raw: std²/mean)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Density Distribution of Raw Focus Scores", fontsize=14)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # Plot 4: Density plot of normalized focus scores with GMM 2D classification
    ax = axes[1, 1]
    # Use same classification as Plot 2 (already determined above)
    ax.hist(
        blurred["focus_score_norm"],
        bins=50,
        alpha=0.7,
        edgecolor="black",
        color="red",
        label="Blurred",
        density=True,
    )
    ax.hist(
        in_focus["focus_score_norm"],
        bins=50,
        alpha=0.7,
        edgecolor="black",
        color="blue",
        label="In-Focus",
        density=True,
    )

    # Add threshold line only if using threshold method
    if not has_gmm_2d:
        ax.axvline(
            x=threshold,
            color="black",
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({threshold})",
        )

    ax.set_xlabel("Focus Score (Normalized)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(
        f"Density Distribution of Normalized Focus Scores (Threshold: {threshold})",
        fontsize=14,
    )
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.savefig(
        figures_dir / "roi_focus_distribution.pdf", dpi=300, bbox_inches="tight"
    )
    plt.savefig(
        figures_dir / "roi_focus_distribution.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV
    df_dist = df_grid_roi[["focus_score", "focus_score_norm"]].copy()
    df_dist["is_low_nuclear_texture"] = df_dist["focus_score_norm"] <= threshold
    df_dist.to_csv(figures_source_dir / "roi_focus_distribution.csv", index=False)

    # Print summary statistics
    logging.info("  Focus score distribution summary:")
    logging.info(f"    Raw focus score - Mean: {mean_fs:.4f}, Median: {median_fs:.4f}")
    logging.info(
        f"    Normalized focus score - Mean: {mean_norm:.4f}, Median: {median_norm:.4f}"
    )
    logging.info(
        f"    Tiles blurred (≤{threshold}): {len(blurred)} ({pct_blurred:.1f}%)"
    )
    logging.info(
        f"    Tiles in-focus (>{threshold}): {len(in_focus)} ({100 - pct_blurred:.1f}%)"
    )


def calculate_roi_intensities(xoa_morphology_files, df_grid_roi):
    """
    Calculate mean intensity per tile for Boundary and IntRNA channels.

    Note: If intensities are already calculated in calculate_roi_focusscore() (new behavior),
    this function will detect that and return the DataFrame unchanged.

    Parameters:
    -----------
    xoa_morphology_files : list
        List of paths to morphology image files
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: roi_id, x1, x2, y1, y2, raw_intensity (DAPI)
        If dapi_intensity, boundary_intensity, intrna_intensity already exist, returns unchanged.

    Returns:
    --------
    pandas DataFrame
        DataFrame with columns:
        - dapi_intensity: Mean DAPI intensity per tile
        - boundary_intensity: Mean boundary intensity per tile (if available)
        - intrna_intensity: Mean IntRNA intensity per tile (if available)
    """

    # Check if intensities are already calculated (new behavior in calculate_roi_focusscore)
    # Note: boundary_intensity may be NaN if only DAPI channel is available, but column should exist
    if "dapi_intensity" in df_grid_roi.columns:
        # Intensities already calculated, just return the DataFrame
        logging.info(
            "  Intensities already calculated in calculate_roi_focusscore(), skipping recalculation"
        )
        return df_grid_roi.copy()

    # Legacy behavior: calculate intensities if not already present
    # Load channels robustly from multi-channel stack or split channel files.
    dapi_image, boundary_image, intrna_image = _load_morphology_channels(
        xoa_morphology_files, level=0
    )
    has_boundary = boundary_image is not None
    has_intrna = intrna_image is not None

    # Create output DataFrame
    df_roi_intensities = df_grid_roi.copy()
    if "dapi_intensity" not in df_roi_intensities.columns:
        df_roi_intensities["dapi_intensity"] = df_roi_intensities.get(
            "raw_intensity", np.nan
        )  # Rename for consistency

    # Calculate intensities for each ROI
    n_rois = len(df_grid_roi)

    # Pre-extract ROI coordinate arrays for vectorized access
    x1_arr = df_grid_roi["x1"].values.astype(int)
    x2_arr = df_grid_roi["x2"].values.astype(int)
    y1_arr = df_grid_roi["y1"].values.astype(int)
    y2_arr = df_grid_roi["y2"].values.astype(int)

    # Boundary intensity
    if has_boundary:
        boundary_intensities = np.empty(n_rois, dtype=np.float64)
        for i in range(n_rois):
            boundary_intensities[i] = np.mean(
                boundary_image[y1_arr[i] : y2_arr[i], x1_arr[i] : x2_arr[i]]
            )
        df_roi_intensities["boundary_intensity"] = boundary_intensities
        del boundary_image
    else:
        df_roi_intensities["boundary_intensity"] = np.nan
        logging.warning(
            "Warning: Boundary channel not available, setting boundary_intensity to NaN"
        )

    # IntRNA intensity
    if has_intrna:
        intrna_intensities = np.empty(n_rois, dtype=np.float64)
        for i in range(n_rois):
            intrna_intensities[i] = np.mean(
                intrna_image[y1_arr[i] : y2_arr[i], x1_arr[i] : x2_arr[i]]
            )
        df_roi_intensities["intrna_intensity"] = intrna_intensities
    else:
        df_roi_intensities["intrna_intensity"] = np.nan
        logging.warning(
            "Warning: IntRNA channel not available, setting intrna_intensity to NaN"
        )

    # Clean up
    del dapi_image

    return df_roi_intensities


def assess_raw_intensity_quality(
    df_roi_intensities,
    dapi_threshold_critical=500,
    boundary_threshold_critical=100,
    intrna_threshold_critical=300,
    min_tissue_coverage=ROI_MIN_TISSUE_COVERAGE_FOR_INTENSITY_QC,
    channel_pct_thresholds=None,
):
    """
    Assess raw intensity quality from per-tile intensities.
    Handles missing channels gracefully (NaN values).

    Note: This function can apply an additional tissue coverage filter before
    intensity QC. If ``tissue_coverage`` is present, only tiles with
    ``tissue_coverage >= min_tissue_coverage`` are used for threshold-based
    percentage calculations. This reduces false FAILs from low-content/background
    edge tiles that technically overlap tissue but are mostly non-tissue.

    Threshold Method:
    -----------------
    Each channel has a single absolute minimum intensity threshold (``intensity_critical``
    from YAML).  The metric is the fraction of tissue tiles whose mean intensity falls
    below that threshold: ``pct_tissue_roi_below_critical``.

    Critical thresholds (defaults):
      - DAPI: 500  — ~0.76% of 16-bit max; typical well-stained range 2000-8000
      - Boundary: 100  — membrane markers 5-10× lower than DAPI
      - IntRNA: 300  — rRNA markers 2-3× lower than DAPI

    Quality Status Determination:
    -----------------------------
    Per-channel verdict is driven by ``pct_tissue_roi_below_critical`` compared
    against YAML ``intensity_warn`` / ``intensity_fail`` fractions (converted to %):
      - **critical**: pct_tissue_roi_below_critical > fail%  OR  mean < critical_threshold
      - **warning**:  pct_tissue_roi_below_critical > warn%
      - **good**: otherwise

    Parameters:
    -----------
    df_roi_intensities : pandas DataFrame
        DataFrame with per-tile intensities (dapi_intensity, boundary_intensity, intrna_intensity).
    dapi_threshold_critical : float, optional
        Critical DAPI threshold (default: 500).
    boundary_threshold_critical : float, optional
        Critical boundary threshold (default: 100).
    intrna_threshold_critical : float, optional
        Critical IntRNA threshold (default: 300).

    Returns:
    --------
    dict
        Per-channel dict with:
        - mean, median, p10, p25, p75, p90: Intensity statistics
        - critical_threshold: Absolute minimum intensity threshold
        - n_tissue_rois_below_critical: Count of tissue tiles below threshold
        - pct_tissue_roi_below_critical: Percentage of tissue tiles below threshold
        - pct_warn_threshold, pct_fail_threshold: YAML-derived population thresholds (%)
        - quality_status: 'good', 'warning', 'critical', or 'not_available'
    """
    stats = {}
    _cpct = channel_pct_thresholds or {}
    df_work = df_roi_intensities
    n_rois_input = int(len(df_roi_intensities))
    if "tissue_coverage" in df_roi_intensities.columns:
        df_work = df_roi_intensities[
            df_roi_intensities["tissue_coverage"] >= float(min_tissue_coverage)
        ]
    n_rois_used = int(len(df_work))
    stats["_intensity_qc_scope"] = {
        "n_rois_input": n_rois_input,
        "n_rois_used": n_rois_used,
        "min_tissue_coverage": float(min_tissue_coverage),
        "uses_tissue_coverage_filter": "tissue_coverage" in df_roi_intensities.columns,
    }

    # YAML key → function-local name mapping
    _yaml_keys = {"dapi": "DAPI", "boundary": "boundary", "intrna": "intRNA"}

    for channel, critical_threshold in [
        ("dapi", dapi_threshold_critical),
        ("boundary", boundary_threshold_critical),
        ("intrna", intrna_threshold_critical),
    ]:
        # Per-channel warn/fail fractions from YAML (or defaults).
        _ch = _cpct.get(_yaml_keys.get(channel, channel)) or {}
        pct_warn_frac = float(_ch.get("intensity_warn", 0.15))
        pct_fail_frac = float(_ch.get("intensity_fail", 0.35))
        col_name = f"{channel}_intensity"
        intensities = df_work[col_name].values

        # Check if channel is available (not all NaN)
        if np.all(np.isnan(intensities)):
            stats[channel] = {
                "mean": np.nan,
                "median": np.nan,
                "p10": np.nan,
                "p25": np.nan,
                "p75": np.nan,
                "p90": np.nan,
                "critical_threshold": float(critical_threshold),
                "n_tissue_rois_below_critical": 0,
                "pct_tissue_roi_below_critical": 0.0,
                "quality_status": "not_available",
            }
            continue

        # Filter out NaN values for calculations
        intensities_valid = intensities[~np.isnan(intensities)]

        if len(intensities_valid) == 0:
            stats[channel] = {
                "mean": np.nan,
                "median": np.nan,
                "p10": np.nan,
                "p25": np.nan,
                "p75": np.nan,
                "p90": np.nan,
                "critical_threshold": float(critical_threshold),
                "n_tissue_rois_below_critical": 0,
                "pct_tissue_roi_below_critical": 0.0,
                "quality_status": "not_available",
            }
            continue

        # Calculate statistics (using valid values only)
        mean_int = np.mean(intensities_valid)
        median_int = np.median(intensities_valid)
        p10 = np.percentile(intensities_valid, 10)
        p25 = np.percentile(intensities_valid, 25)
        p75 = np.percentile(intensities_valid, 75)
        p90 = np.percentile(intensities_valid, 90)

        # Count tissue tiles below the single critical intensity threshold
        n_below = int(np.sum(intensities_valid < critical_threshold))
        pct_tissue_roi_below_critical = (n_below / len(intensities_valid)) * 100

        # Determine quality status (thresholds from YAML, as percentages 0-100)
        pct_fail = pct_fail_frac * 100.0
        pct_warn = pct_warn_frac * 100.0
        if pct_tissue_roi_below_critical > pct_fail or mean_int < critical_threshold:
            quality_status = "critical"
        elif pct_tissue_roi_below_critical > pct_warn:
            quality_status = "warning"
        else:
            quality_status = "good"

        stats[channel] = {
            "mean": float(mean_int),
            "median": float(median_int),
            "p10": float(p10),
            "p25": float(p25),
            "p75": float(p75),
            "p90": float(p90),
            "critical_threshold": float(critical_threshold),
            "pct_warn_threshold": float(pct_warn),
            "pct_fail_threshold": float(pct_fail),
            "n_tissue_rois_below_critical": n_below,
            "pct_tissue_roi_below_critical": float(pct_tissue_roi_below_critical),
            "quality_status": quality_status,
        }

    # Overall quality
    statuses = [stats[ch]["quality_status"] for ch in ["dapi", "boundary", "intrna"]]
    if n_rois_used == 0:
        overall_quality = "not_available"
    elif "critical" in statuses or statuses.count("warning") >= 2:
        overall_quality = "critical"
    elif "warning" in statuses:
        overall_quality = "warning"
    else:
        overall_quality = "good"

    stats["overall_quality"] = overall_quality

    return stats


def plot_focus_score_vs_laplacian(df_grid_roi, figures_dir, figures_source_dir):
    """
    Compare Laplacian variance vs original focus score (std²/mean).

    Creates scatter plot and distribution histograms comparing the two focus metrics.

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: dapi_focus_score, dapi_lap_var
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    """
    # Check if required columns exist
    if "dapi_lap_var" not in df_grid_roi.columns:
        logging.warning(
            "  Warning: dapi_lap_var column not found, skipping Laplacian comparison plots"
        )
        return

    if "dapi_focus_score" not in df_grid_roi.columns:
        # Fallback to generic focus_score
        if "focus_score" not in df_grid_roi.columns:
            logging.warning(
                "  Warning: No focus score column found, skipping Laplacian comparison plots"
            )
            return
        focus_col = "focus_score"
    else:
        focus_col = "dapi_focus_score"

    # Filter out NaN values
    valid_mask = df_grid_roi["dapi_lap_var"].notna() & df_grid_roi[focus_col].notna()
    df_valid = df_grid_roi[valid_mask].copy()

    if len(df_valid) == 0:
        logging.warning(
            "  Warning: No valid data for Laplacian comparison, skipping plots"
        )
        return

    # Scatter plot: log1p(focus_score) vs log1p(laplacian_variance)
    logging.info("Generating Figure: Focus score vs Laplacian variance comparison...")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    # Color by GMM 2D classification if available
    has_gmm_2d = "is_blurred_gmm_2d" in df_valid.columns
    focus_vals = df_valid[focus_col].values
    lap_vals = df_valid["dapi_lap_var"].values
    if has_gmm_2d:
        # Filter to valid mask for GMM 2D
        valid_gmm_mask = df_valid["is_blurred_gmm_2d"].notna().values
        if valid_gmm_mask.sum() > 0:
            sidx = _subsample_idx(valid_gmm_mask)
            is_blurred_sub = df_valid["is_blurred_gmm_2d"].values[sidx].astype(bool)
            colors = np.where(is_blurred_sub, "red", "blue")
            ax.scatter(
                np.log1p(focus_vals[sidx]),
                np.log1p(lap_vals[sidx]),
                s=5,
                alpha=0.3,
                c=colors,
                rasterized=True,
            )
            from matplotlib.patches import Patch

            legend_elements = [
                Patch(facecolor="blue", label="In-Focus (2D GMM)"),
                Patch(facecolor="red", label="Blurred (2D GMM)"),
            ]
            ax.legend(handles=legend_elements, fontsize=10)
        else:
            sidx = _subsample_idx(np.ones(len(df_valid), dtype=bool))
            ax.scatter(
                np.log1p(focus_vals[sidx]),
                np.log1p(lap_vals[sidx]),
                s=5,
                alpha=0.3,
                rasterized=True,
            )
    else:
        sidx = _subsample_idx(np.ones(len(df_valid), dtype=bool))
        ax.scatter(
            np.log1p(focus_vals[sidx]),
            np.log1p(lap_vals[sidx]),
            s=5,
            alpha=0.3,
            rasterized=True,
        )

    ax.set_xlabel("log1p(std²/mean) (DAPI)", fontsize=12)
    ax.set_ylabel("log1p(Laplacian variance) (DAPI)", fontsize=12)
    if has_gmm_2d:
        ax.set_title(
            "Focus Score vs Laplacian Variance (2D GMM Classification)", fontsize=14
        )
    else:
        ax.set_title("Focus Score vs Laplacian Variance", fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    # Save to methodology assessment folder (subfolder of figures_dir)
    methodology_figures_dir = figures_dir / "figures_methodology_assessment"
    methodology_figures_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        methodology_figures_dir / "focus_score_vs_laplacian.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        methodology_figures_dir / "focus_score_vs_laplacian.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Distribution comparison
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].hist(df_valid[focus_col], bins=50, alpha=0.7, edgecolor="black")
    axes[0].set_title("std²/mean (DAPI Focus Score)", fontsize=12)
    axes[0].set_xlabel("Focus Score (std²/mean)", fontsize=11)
    axes[0].set_ylabel("Number of tiles", fontsize=11)
    axes[0].grid(True, alpha=0.3, axis="y")

    axes[1].hist(df_valid["dapi_lap_var"], bins=50, alpha=0.7, edgecolor="black")
    axes[1].set_title("Laplacian Variance (DAPI)", fontsize=12)
    axes[1].set_xlabel("Laplacian Variance", fontsize=11)
    axes[1].set_ylabel("Number of tiles", fontsize=11)
    axes[1].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    # Save to methodology assessment folder
    plt.savefig(
        methodology_figures_dir / "focus_score_vs_laplacian_distributions.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        methodology_figures_dir / "focus_score_vs_laplacian_distributions.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save source data
    df_comparison = df_valid[[focus_col, "dapi_lap_var"]].copy()
    df_comparison["log1p_focus_score"] = np.log1p(df_comparison[focus_col])
    df_comparison["log1p_lap_var"] = np.log1p(df_comparison["dapi_lap_var"])
    df_comparison.to_csv(
        figures_source_dir / "focus_score_vs_laplacian.csv", index=False
    )

    # Calculate correlation
    correlation = np.corrcoef(
        np.log1p(df_valid[focus_col]), np.log1p(df_valid["dapi_lap_var"])
    )[0, 1]
    logging.info(f"  Correlation (log1p): {correlation:.4f}")


def plot_intensity_assessment(
    df_roi_intensities, intensity_stats, small0, figures_dir, figures_source_dir
):
    """
    Create visualization of per-tile intensity distributions and spatial heatmaps.

    Note: Intensity distributions are calculated from tissue-filtered tiles only
    (tiles with any tissue overlap, i.e., coverage > 0%). Background/empty regions
    of the image are excluded from the distributions.

    Parameters:
    -----------
    df_roi_intensities : pandas DataFrame
        DataFrame with per-tile intensities and coordinates (tissue-filtered tiles only)
    intensity_stats : dict
        Dictionary from assess_raw_intensity_quality()
    small0 : numpy.ndarray
        Downsampled DAPI image (for spatial overlay)
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    """
    downsample_factor = 8
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Calculate ROI centers
    df_roi_intensities["center_x"] = (
        df_roi_intensities["x1"] + df_roi_intensities["x2"]
    ) / 2
    df_roi_intensities["center_y"] = (
        df_roi_intensities["y1"] + df_roi_intensities["y2"]
    ) / 2

    channels = [
        ("dapi", "DAPI", intensity_stats["dapi"]),
        ("boundary", "Boundary", intensity_stats["boundary"]),
        ("intrna", "IntRNA", intensity_stats["intrna"]),
    ]

    for col_idx, (channel, channel_name, stats) in enumerate(channels):
        intensity_col = f"{channel}_intensity"
        intensities = df_roi_intensities[intensity_col].values

        # Check if channel is available
        is_available = not np.all(np.isnan(intensities))
        intensities_valid = (
            intensities[~np.isnan(intensities)] if is_available else np.array([])
        )

        # Row 1: Distribution plots
        ax = axes[0, col_idx]

        if is_available and len(intensities_valid) > 0:
            ax.hist(intensities_valid, bins=50, alpha=0.7, edgecolor="black")

            # Add threshold lines (only if thresholds are valid)
            if not np.isnan(stats["critical_threshold"]):
                ax.axvline(
                    stats["critical_threshold"],
                    color="red",
                    linestyle="--",
                    linewidth=2,
                    label=f"Min intensity ({stats['critical_threshold']:.0f})",
                )

            # No green "optimal range" band — only the red critical threshold
            # line (from YAML intensity_critical) is shown to avoid implying an
            # uncalibrated optimal window.

            # Add statistics text
            if not np.isnan(stats["mean"]):
                stats_text = f"Mean: {stats['mean']:.0f}\nMedian: {stats['median']:.0f}\nP10: {stats['p10']:.0f}\nP90: {stats['p90']:.0f}"
                ax.text(
                    0.98,
                    0.98,
                    stats_text,
                    transform=ax.transAxes,
                    verticalalignment="top",
                    horizontalalignment="right",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
                    fontsize=9,
                )
        else:
            ax.text(
                0.5,
                0.5,
                f"{channel_name} channel\nnot available",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=14,
                bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
            )

        ax.set_xlabel(f"{channel_name} Intensity", fontsize=12)
        ax.set_ylabel("Number of tiles", fontsize=12)
        status_text = (
            stats["quality_status"].upper()
            if stats["quality_status"] != "not_available"
            else "NOT AVAILABLE"
        )
        ax.set_title(
            f"{channel_name} Intensity Distribution\n(Status: {status_text})",
            fontsize=14,
        )
        if is_available and len(intensities_valid) > 0:
            ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Row 2: Spatial heatmaps
        ax = axes[1, col_idx]

        # Get image dimensions for setting axes limits
        img_height, img_width = small0.shape

        # Show downsampled DAPI as background (optional, light)
        ax.imshow(
            small0,
            cmap="Greys_r",
            vmax=np.percentile(small0, 99),
            alpha=0.3,
            aspect="auto",
            origin="upper",
            extent=[0, img_width, img_height, 0],
        )

        # Scale coordinates to downsampled space
        x_scaled = df_roi_intensities["center_x"] / downsample_factor
        y_scaled = df_roi_intensities["center_y"] / downsample_factor

        if is_available and len(intensities_valid) > 0:
            # Create scatter plot heatmap (only for valid intensities)
            valid_mask = ~np.isnan(intensities)

            # Filter coordinates to be within image bounds
            x_vals = x_scaled.values if hasattr(x_scaled, "values") else x_scaled
            y_vals = y_scaled.values if hasattr(y_scaled, "values") else y_scaled
            in_bounds = (
                (x_vals >= 0)
                & (x_vals < img_width)
                & (y_vals >= 0)
                & (y_vals < img_height)
            )
            valid_mask = valid_mask & in_bounds

            if valid_mask.sum() > 0:
                # Use hexbin for spatial heatmap — bins ALL data, renders O(bins) not O(N)
                hb = ax.hexbin(
                    x_vals[valid_mask],
                    y_vals[valid_mask],
                    C=intensities[valid_mask],
                    reduce_C_function=np.mean,
                    gridsize=100,
                    cmap="viridis",
                    mincnt=1,
                )

                # Add colorbar
                cbar = plt.colorbar(hb, ax=ax)
                cbar.set_label(f"{channel_name} Intensity", fontsize=10)
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No valid data points\nwithin image bounds",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=12,
                    bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
                )
        else:
            ax.text(
                0.5,
                0.5,
                f"{channel_name} channel\nnot available",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=14,
                bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
            )

        # Set axes limits explicitly
        ax.set_xlim(0, img_width)
        ax.set_ylim(img_height, 0)  # Reversed because origin='upper'
        ax.set_xlabel("X coordinate (downsampled)", fontsize=12)
        ax.set_ylabel("Y coordinate (downsampled)", fontsize=12)
        ax.set_title(f"{channel_name} Intensity Spatial Heatmap", fontsize=14)
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(figures_dir / "intensity_assessment.pdf", dpi=300, bbox_inches="tight")
    plt.savefig(figures_dir / "intensity_assessment.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save source data
    df_roi_intensities.to_csv(
        figures_source_dir / "intensity_assessment.csv", index=False
    )

    # Save statistics summary
    df_stats = pd.DataFrame(
        {
            "channel": ["dapi", "boundary", "intrna"],
            "mean_intensity": [
                intensity_stats["dapi"]["mean"],
                intensity_stats["boundary"]["mean"],
                intensity_stats["intrna"]["mean"],
            ],
            "median_intensity": [
                intensity_stats["dapi"]["median"],
                intensity_stats["boundary"]["median"],
                intensity_stats["intrna"]["median"],
            ],
            "p10_intensity": [
                intensity_stats["dapi"]["p10"],
                intensity_stats["boundary"]["p10"],
                intensity_stats["intrna"]["p10"],
            ],
            "p90_intensity": [
                intensity_stats["dapi"]["p90"],
                intensity_stats["boundary"]["p90"],
                intensity_stats["intrna"]["p90"],
            ],
            "critical_threshold": [
                intensity_stats["dapi"]["critical_threshold"],
                intensity_stats["boundary"]["critical_threshold"],
                intensity_stats["intrna"]["critical_threshold"],
            ],
            "pct_tissue_roi_below_critical": [
                intensity_stats["dapi"]["pct_tissue_roi_below_critical"],
                intensity_stats["boundary"]["pct_tissue_roi_below_critical"],
                intensity_stats["intrna"]["pct_tissue_roi_below_critical"],
            ],
            "quality_status": [
                intensity_stats["dapi"]["quality_status"],
                intensity_stats["boundary"]["quality_status"],
                intensity_stats["intrna"]["quality_status"],
            ],
        }
    )
    df_stats.to_csv(
        figures_source_dir / "intensity_assessment_statistics.csv", index=False
    )


def generate_roi_figures(
    data,
    small0,
    small1,
    small2,
    distance_map,
    distance_map2,
    whole_sample,
    holes,
    dense_intensity_regions,
    df_grid_roi,
    df_roi_intensities,
    intensity_stats,
    xoa_morphology_files,
    focus_maps=None,
    snr_thresholds=None,
    blur_prob_threshold=0.5,
):
    """
    Generate cell-independent tile-based figures using multithreading.

    Parameters:
    -----------
    data : dict
        Dictionary with 'figures_dir' key
    small0, small1, small2 : numpy.ndarray
        Downsampled morphology images
    distance_map, distance_map2 : numpy.ndarray
        Distance maps (edge and holes)
    whole_sample, holes, dense_intensity_regions : numpy.ndarray
        Tissue masks
    df_grid_roi : pandas DataFrame
        Grid ROI focus scores
    df_roi_intensities : pandas DataFrame
        Tile intensity measurements
    intensity_stats : dict
        Intensity quality assessment statistics
    xoa_morphology_files : list
        List of morphology file paths
    focus_maps : dict, optional
        Focus map arrays for heatmap overlay
    snr_thresholds : dict, optional
        YAML ``snr.roi_tx`` block for SNR heatmap threshold lines
    """
    figures_dir = data["figures_dir"]
    figures_source_dir = figures_dir / "figures_source"
    figures_source_dir.mkdir(parents=True, exist_ok=True)
    # Create methodology assessment folder for comparison figures
    methodology_figures_dir = figures_dir / "figures_methodology_assessment"
    methodology_figures_dir.mkdir(parents=True, exist_ok=True)

    def _fig1_distance_edge():
        logging.info("Generating Figure 1: Distance map (edge)...")
        t_fig = time.time()
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(distance_map, cmap="gray")
        ax.set_title("Distance to Edge")
        plt.tight_layout()
        plt.savefig(figures_dir / "distance_map_edge.pdf", dpi=300, bbox_inches="tight")
        plt.savefig(figures_dir / "distance_map_edge.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(distance_map).to_csv(
            figures_source_dir / "distance_map_edge.csv", index=False, header=False
        )
        logging.info(
            f"[TIMING] Figure 1 (distance map edge): {time.time() - t_fig:.1f}s"
        )

    def _fig2_distance_holes():
        logging.info("Generating Figure 2: Distance map (holes)...")
        t_fig = time.time()
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(distance_map2, cmap="gray")
        ax.set_title("Distance to Nearest Hole")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "distance_map_holes.pdf", dpi=300, bbox_inches="tight"
        )
        plt.savefig(
            figures_dir / "distance_map_holes.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        pd.DataFrame(distance_map2).to_csv(
            figures_source_dir / "distance_map_holes.csv", index=False, header=False
        )
        logging.info(
            f"[TIMING] Figure 2 (distance map holes): {time.time() - t_fig:.1f}s"
        )

    def _fig3_morphology_overview():
        logging.info("Generating Figure 3: Morphology overview...")
        t_fig = time.time()
        r = (small0.shape[0] / small0.shape[1]) * 12
        fig, ax = plt.subplots(2, 2, figsize=(12, r))
        ax[0, 0].imshow(small0, cmap="Greys_r", vmax=np.percentile(small0, 99))
        ax[0, 0].set_title("DAPI")
        ax[0, 0].set_aspect("equal")
        if small1 is not None:
            ax[0, 1].imshow(small1, cmap="Greys_r", vmax=np.percentile(small1, 99))
            ax[0, 1].set_title("Boundary")
            ax[0, 1].set_aspect("equal")
        else:
            ax[0, 1].set_title("Boundary (not available)")
            ax[0, 1].axis("off")
        if small2 is not None:
            ax[1, 0].imshow(small2, cmap="Greys_r", vmax=np.percentile(small2, 99))
            ax[1, 0].set_title("Interior")
            ax[1, 0].set_aspect("equal")
        else:
            ax[1, 0].set_title("Interior (not available)")
            ax[1, 0].axis("off")
        ax[1, 1].set_title("Sample dense intensity regions map")
        ax[1, 1].imshow(color.label2rgb(dense_intensity_regions, bg_label=0))
        ax[1, 1].set_aspect("equal")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "morphology_overview.pdf", dpi=300, bbox_inches="tight"
        )
        plt.savefig(
            figures_dir / "morphology_overview.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        pd.DataFrame(small0).to_csv(
            figures_source_dir / "morphology_overview_DAPI.csv",
            index=False,
            header=False,
        )
        if small1 is not None:
            pd.DataFrame(small1).to_csv(
                figures_source_dir / "morphology_overview_Boundary.csv",
                index=False,
                header=False,
            )
        if small2 is not None:
            pd.DataFrame(small2).to_csv(
                figures_source_dir / "morphology_overview_Interior.csv",
                index=False,
                header=False,
            )
        pd.DataFrame(dense_intensity_regions).to_csv(
            figures_source_dir / "morphology_overview_DenseIntensityRegions.csv",
            index=False,
            header=False,
        )
        logging.info(
            f"[TIMING] Figure 3 (morphology overview): {time.time() - t_fig:.1f}s"
        )

    def _fig4_imageqc_masks():
        logging.info("Generating Figure 4: ImageQC masks...")
        t_fig = time.time()
        fig, ax = plt.subplots(2, 2, figsize=(12, 10))
        ax[0, 0].set_title("DAPI morphology image")
        ax[0, 0].imshow(small0, cmap="Greys_r", vmax=np.percentile(small0, 99))
        ax[0, 0].axis("off")
        ax[0, 1].set_title("Whole Sample mask")
        ax[0, 1].imshow(color.label2rgb(whole_sample, bg_label=0))
        ax[0, 1].axis("off")
        ax[1, 0].set_title("Holes in sample")
        ax[1, 0].imshow(color.label2rgb(holes, bg_label=0))
        ax[1, 0].axis("off")
        ax[1, 1].set_title("Sample dense intensity regions")
        ax[1, 1].imshow(color.label2rgb(dense_intensity_regions, bg_label=0))
        ax[1, 1].axis("off")
        plt.tight_layout()
        plt.savefig(figures_dir / "imageqc_masks.pdf", dpi=300, bbox_inches="tight")
        plt.savefig(figures_dir / "imageqc_masks.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(small0).to_csv(
            figures_source_dir / "imageqc_masks_DAPI.csv", index=False, header=False
        )
        pd.DataFrame(whole_sample).to_csv(
            figures_source_dir / "imageqc_masks_WholeSample.csv",
            index=False,
            header=False,
        )
        pd.DataFrame(holes).to_csv(
            figures_source_dir / "imageqc_masks_Holes.csv", index=False, header=False
        )
        pd.DataFrame(dense_intensity_regions).to_csv(
            figures_source_dir / "imageqc_masks_DenseIntensityRegions.csv",
            index=False,
            header=False,
        )
        logging.info(f"[TIMING] Figure 4 (ImageQC masks): {time.time() - t_fig:.1f}s")

    def _fig5_focus_heatmap():
        logging.info("Generating Figure 5: Grid tile focus heatmap...")
        t_fig = time.time()
        plot_grid_roi_focus_heatmap(
            df_grid_roi,
            small0,
            figures_dir,
            figures_source_dir,
            threshold=-1.0,
            focus_maps=focus_maps,
        )
        logging.info(f"[TIMING] Figure 5 (focus heatmap): {time.time() - t_fig:.1f}s")

    def _fig5a_blur_heatmap():
        prob_col = (
            "blur_prob_gmm_2d"
            if "blur_prob_gmm_2d" in df_grid_roi.columns
            else "blur_prob_gmm"
            if "blur_prob_gmm" in df_grid_roi.columns
            else None
        )
        class_col = (
            "is_blurred_gmm_2d"
            if "is_blurred_gmm_2d" in df_grid_roi.columns
            else "is_blurred_gmm"
            if "is_blurred_gmm" in df_grid_roi.columns
            else None
        )
        if prob_col is None and class_col is None:
            logging.info("Skipping blur heatmap: no blur columns in df_grid_roi")
            return
        logging.info("Generating Figure 5a: Spatial blur heatmap...")
        t_fig = time.time()
        plot_grid_roi_blur_heatmap(
            df_grid_roi,
            small0,
            figures_dir,
            figures_source_dir,
            blur_prob_threshold=blur_prob_threshold,
        )
        logging.info(f"[TIMING] Figure 5a (blur heatmap): {time.time() - t_fig:.1f}s")

    def _fig5b_focus_vs_intensity():
        logging.info("Generating Figure 5b: Tile focus score vs intensity...")
        t_fig = time.time()
        plot_roi_focus_vs_intensity(
            df_grid_roi, figures_dir, figures_source_dir, threshold=-1.0
        )
        logging.info(
            f"[TIMING] Figure 5b (focus vs intensity): {time.time() - t_fig:.1f}s"
        )

    def _fig5c_focus_distribution():
        logging.info("Generating Figure 5c: Tile focus score distribution...")
        t_fig = time.time()
        plot_roi_focus_distribution(
            df_grid_roi, figures_dir, figures_source_dir, threshold=-1.0
        )
        logging.info(
            f"[TIMING] Figure 5c (focus distribution): {time.time() - t_fig:.1f}s"
        )

    def _fig5d_focus_vs_laplacian():
        logging.info(
            "Generating Figure 5d: Focus score vs Laplacian variance comparison..."
        )
        t_fig = time.time()
        plot_focus_score_vs_laplacian(df_grid_roi, figures_dir, figures_source_dir)
        logging.info(
            f"[TIMING] Figure 5d (focus vs Laplacian): {time.time() - t_fig:.1f}s"
        )

    def _fig6_intensity_assessment():
        logging.info("Generating Figure 6: Intensity assessment...")
        t_fig = time.time()
        plot_intensity_assessment(
            df_roi_intensities, intensity_stats, small0, figures_dir, figures_source_dir
        )
        logging.info(
            f"[TIMING] Figure 6 (intensity assessment): {time.time() - t_fig:.1f}s"
        )

    def _fig7_snr_heatmap():
        if "neg_pct" not in df_grid_roi.columns:
            logging.info("Skipping SNR heatmap: no SNR columns in df_grid_roi")
            return
        logging.info("Generating Figure 7: SNR spatial heatmap...")
        t_fig = time.time()
        plot_snr_roi_heatmap(
            df_grid_roi,
            small0,
            figures_dir,
            figures_source_dir,
        )
        logging.info(f"[TIMING] Figure 7 (SNR heatmap): {time.time() - t_fig:.1f}s")

    def _fig8_concordance():
        if "roi_tx_snr_ratio" not in df_grid_roi.columns:
            logging.info("Skipping concordance plot: no SNR columns in df_grid_roi")
            return
        logging.info("Generating Figure 8: Cross-section concordance...")
        t_fig = time.time()
        plot_cross_section_concordance(
            df_grid_roi,
            figures_dir,
            figures_source_dir,
            snr_thresholds=snr_thresholds,
        )
        logging.info(f"[TIMING] Figure 8 (concordance): {time.time() - t_fig:.1f}s")

    tasks = [
        _fig1_distance_edge,
        _fig2_distance_holes,
        _fig3_morphology_overview,
        _fig4_imageqc_masks,
        _fig5_focus_heatmap,
        _fig5a_blur_heatmap,
        _fig5b_focus_vs_intensity,
        _fig5c_focus_distribution,
        _fig5d_focus_vs_laplacian,
        _fig6_intensity_assessment,
        _fig7_snr_heatmap,
        _fig8_concordance,
    ]

    # Use multiprocessing.Process with fork context -- matplotlib is NOT thread-safe.
    # fork is intentional: this code runs exclusively in Linux/Docker containers where
    # fork is safe, and the closures capture large numpy arrays that cannot be pickled
    # (required by spawn/forkserver). Do NOT change to spawn/forkserver without
    # restructuring the closure-based task definitions.
    max_workers = min(len(tasks), os.cpu_count() or 4)
    mp_ctx = multiprocessing.get_context("fork")

    # Run in batches to limit concurrent processes
    failed = []
    for batch_start in range(0, len(tasks), max_workers):
        batch = tasks[batch_start : batch_start + max_workers]
        processes = []
        for fn in batch:
            p = mp_ctx.Process(target=_run_figure_task, args=(fn,))
            p.start()
            processes.append((p, fn.__name__))
        for p, name in processes:
            p.join()
            if p.exitcode != 0:
                failed.append(name)

    if failed:
        raise RuntimeError(f"ROI figure generation failed for: {', '.join(failed)}")

    logging.info("All tile-based figures generated successfully!")
    logging.info(f"Saved figures to {figures_dir}")
    logging.info(f"Saved source data to {figures_source_dir}")


def save_roi_qc_metrics(
    df_grid_roi,
    intensity_stats,
    outdir,
    roi_size=None,
    snr_summary=None,
    distance_map=None,
    distance_map2=None,
    edge_distance_threshold: float = -25.0,
    hole_distance_threshold: float = -25.0,
    min_tissue_coverage_for_qc: float = ROI_MIN_TISSUE_COVERAGE_FOR_INTENSITY_QC,
    qc_thresholds: dict | None = None,
    lap_sigma: float | None = None,
):
    """
    Save tile-based QC metrics to JSON file.

    Parameters
    ----------
    snr_summary : dict or None
        If provided, stored under ``snr`` (from :mod:`snr_metrics`).
    """
    import json

    # Existing stats
    roi_metrics = {
        "roi_size_pixels": roi_size,
        "total_rois": len(df_grid_roi),
        "rois_in_tissue": int((df_grid_roi["tissue_coverage"] > 0).sum()),
        "focus_score": {
            "mean": float(df_grid_roi["focus_score"].mean()),
            "median": float(df_grid_roi["focus_score"].median()),
            "std": float(df_grid_roi["focus_score"].std()),
            "min": float(df_grid_roi["focus_score"].min()),
            "max": float(df_grid_roi["focus_score"].max()),
            "tissue_median": float(
                df_grid_roi.loc[
                    df_grid_roi["tissue_coverage"] >= min_tissue_coverage_for_qc,
                    "focus_score",
                ].median()
            )
            if "tissue_coverage" in df_grid_roi.columns
            else float(df_grid_roi["focus_score"].median()),
        },
        "focus_score_norm": {
            "mean": float(df_grid_roi["focus_score_norm"].mean()),
            "median": float(df_grid_roi["focus_score_norm"].median()),
            "std": float(df_grid_roi["focus_score_norm"].std()),
            "min": float(df_grid_roi["focus_score_norm"].min()),
            "max": float(df_grid_roi["focus_score_norm"].max()),
        },
        "raw_intensity": {
            "mean": float(df_grid_roi["raw_intensity"].mean()),
            "median": float(df_grid_roi["raw_intensity"].median()),
            "std": float(df_grid_roi["raw_intensity"].std()),
            "min": float(df_grid_roi["raw_intensity"].min()),
            "max": float(df_grid_roi["raw_intensity"].max()),
        },
        "tissue_coverage": {
            "mean": float(df_grid_roi["tissue_coverage"].mean()),
            "median": float(df_grid_roi["tissue_coverage"].median()),
            "min": float(df_grid_roi["tissue_coverage"].min()),
            "max": float(df_grid_roi["tissue_coverage"].max()),
        },
        "intensity_quality": intensity_stats,
    }
    total_rois = int(len(df_grid_roi))
    rois_in_tissue = int((df_grid_roi["tissue_coverage"] > 0).sum())
    roi_metrics["tissue_mask_qc"] = {
        "tissue_mask_generated": bool(rois_in_tissue > 0),
        "status": "PASS" if rois_in_tissue > 0 else "FAIL",
        "rois_in_tissue": rois_in_tissue,
        "total_rois": total_rois,
        "tissue_roi_fraction": float(rois_in_tissue / total_rois)
        if total_rois > 0
        else 0.0,
    }

    # Optional: add GMM-based blur summary if columns are present
    if (
        "is_blurred_gmm" in df_grid_roi.columns
        and "is_low_intensity" in df_grid_roi.columns
    ):
        total_rois = len(df_grid_roi)
        n_blurred = int(df_grid_roi["is_blurred_gmm"].sum())
        n_low_int = int(df_grid_roi["is_low_intensity"].sum())
        if "tissue_coverage" in df_grid_roi.columns:
            tissue_mask_qc = (
                df_grid_roi["tissue_coverage"] >= min_tissue_coverage_for_qc
            )
            total_rois_tissue = int(tissue_mask_qc.sum())
            n_blurred_tissue = int(
                df_grid_roi.loc[tissue_mask_qc, "is_blurred_gmm"].sum()
            )
            n_low_int_tissue = int(
                df_grid_roi.loc[tissue_mask_qc, "is_low_intensity"].sum()
            )
        else:
            total_rois_tissue = total_rois
            n_blurred_tissue = n_blurred
            n_low_int_tissue = n_low_int

        roi_metrics["blur_gmm_1d"] = {
            "total_rois": total_rois,
            "rois_blurred_gmm": n_blurred,
            "rois_low_intensity": n_low_int,
            "pct_blurred_gmm": float(n_blurred / total_rois * 100.0)
            if total_rois > 0
            else 0.0,
            "pct_low_intensity": float(n_low_int / total_rois * 100.0)
            if total_rois > 0
            else 0.0,
            # Tissue-aware denominator for report-level interpretation
            "total_rois_tissue_filtered": total_rois_tissue,
            "rois_blurred_gmm_tissue_filtered": n_blurred_tissue,
            "rois_low_intensity_tissue_filtered": n_low_int_tissue,
            "pct_blurred_gmm_tissue_filtered": float(
                n_blurred_tissue / total_rois_tissue * 100.0
            )
            if total_rois_tissue > 0
            else 0.0,
            "pct_low_intensity_tissue_filtered": float(
                n_low_int_tissue / total_rois_tissue * 100.0
            )
            if total_rois_tissue > 0
            else 0.0,
            "tissue_coverage_min_for_qc": float(min_tissue_coverage_for_qc),
        }

        if "blur_prob_gmm" in df_grid_roi.columns:
            valid_probs = df_grid_roi["blur_prob_gmm"].dropna()
            if len(valid_probs) > 0:
                roi_metrics["blur_gmm_1d"]["blur_prob_mean"] = float(valid_probs.mean())
                roi_metrics["blur_gmm_1d"]["blur_prob_median"] = float(
                    valid_probs.median()
                )

    # Optional: add 2D GMM-based blur summary if columns are present
    if "is_blurred_gmm_2d" in df_grid_roi.columns:
        total_rois = len(df_grid_roi)
        n_blurred_2d = int(df_grid_roi["is_blurred_gmm_2d"].sum())
        n_low_int = (
            int(df_grid_roi["is_low_intensity"].sum())
            if "is_low_intensity" in df_grid_roi.columns
            else 0
        )
        if "tissue_coverage" in df_grid_roi.columns:
            tissue_mask_qc = (
                df_grid_roi["tissue_coverage"] >= min_tissue_coverage_for_qc
            )
            total_rois_tissue = int(tissue_mask_qc.sum())
            n_blurred_tissue = int(
                df_grid_roi.loc[tissue_mask_qc, "is_blurred_gmm_2d"].sum()
            )
            n_low_int_tissue = (
                int(df_grid_roi.loc[tissue_mask_qc, "is_low_intensity"].sum())
                if "is_low_intensity" in df_grid_roi.columns
                else 0
            )
        else:
            total_rois_tissue = total_rois
            n_blurred_tissue = n_blurred_2d
            n_low_int_tissue = n_low_int

        roi_metrics["blur_gmm_2d"] = {
            "total_rois": total_rois,
            "rois_blurred_gmm": n_blurred_2d,
            "rois_low_intensity": n_low_int,
            "pct_blurred_gmm": float(n_blurred_2d / total_rois * 100.0)
            if total_rois > 0
            else 0.0,
            "pct_low_intensity": float(n_low_int / total_rois * 100.0)
            if total_rois > 0
            else 0.0,
            # Tissue-aware denominator for report-level interpretation
            "total_rois_tissue_filtered": total_rois_tissue,
            "rois_blurred_gmm_tissue_filtered": n_blurred_tissue,
            "rois_low_intensity_tissue_filtered": n_low_int_tissue,
            "pct_blurred_gmm_tissue_filtered": float(
                n_blurred_tissue / total_rois_tissue * 100.0
            )
            if total_rois_tissue > 0
            else 0.0,
            "pct_low_intensity_tissue_filtered": float(
                n_low_int_tissue / total_rois_tissue * 100.0
            )
            if total_rois_tissue > 0
            else 0.0,
            "tissue_coverage_min_for_qc": float(min_tissue_coverage_for_qc),
        }

        if "blur_prob_gmm_2d" in df_grid_roi.columns:
            valid_probs = df_grid_roi["blur_prob_gmm_2d"].dropna()
            if len(valid_probs) > 0:
                roi_metrics["blur_gmm_2d"]["blur_prob_mean"] = float(valid_probs.mean())
                roi_metrics["blur_gmm_2d"]["blur_prob_median"] = float(
                    valid_probs.median()
                )

        # Compare 1D vs 2D if both exist
        if "is_blurred_gmm" in df_grid_roi.columns:
            both_blurred = (
                df_grid_roi["is_blurred_gmm"] & df_grid_roi["is_blurred_gmm_2d"]
            ).sum()
            both_in_focus = (
                (~df_grid_roi["is_blurred_gmm"]) & (~df_grid_roi["is_blurred_gmm_2d"])
            ).sum()
            agreement = (
                (both_blurred + both_in_focus) / total_rois * 100.0
                if total_rois > 0
                else 0.0
            )
            roi_metrics["gmm_comparison"] = {
                "agreement_pct": float(agreement),
                "both_blurred": int(both_blurred),
                "both_in_focus": int(both_in_focus),
                "only_1d_blurred": int(
                    (
                        df_grid_roi["is_blurred_gmm"]
                        & ~df_grid_roi["is_blurred_gmm_2d"]
                    ).sum()
                ),
                "only_2d_blurred": int(
                    (
                        ~df_grid_roi["is_blurred_gmm"]
                        & df_grid_roi["is_blurred_gmm_2d"]
                    ).sum()
                ),
            }

    if snr_summary is not None:
        roi_metrics["snr"] = snr_summary

    # Optional: morphology ROI summary metrics (vectorized, memory-friendly).
    if (
        distance_map is not None
        and distance_map2 is not None
        and {"x1", "x2", "y1", "y2", "tissue_coverage"}.issubset(df_grid_roi.columns)
    ):
        x1 = df_grid_roi["x1"].to_numpy(np.int64, copy=False)
        x2 = df_grid_roi["x2"].to_numpy(np.int64, copy=False)
        y1 = df_grid_roi["y1"].to_numpy(np.int64, copy=False)
        y2 = df_grid_roi["y2"].to_numpy(np.int64, copy=False)
        tissue_cov = df_grid_roi["tissue_coverage"].to_numpy(np.float64, copy=False)

        # Sample distance maps at ROI centroids in downsampled coordinates.
        cx_ds = np.clip(((x1 + x2) // 2) // 8, 0, distance_map.shape[1] - 1)
        cy_ds = np.clip(((y1 + y2) // 2) // 8, 0, distance_map.shape[0] - 1)
        dist_edge = distance_map[cy_ds, cx_ds]
        dist_hole = distance_map2[cy_ds, cx_ds]

        tissue_mask = tissue_cov > 0.0
        qc_tissue_mask = tissue_cov >= float(min_tissue_coverage_for_qc)
        n_tissue = int(tissue_mask.sum())
        n_qc_tissue = int(qc_tissue_mask.sum())

        edge_zone_frac = (
            float(
                (tissue_mask & (dist_edge > float(edge_distance_threshold))).sum()
                / n_tissue
            )
            if n_tissue > 0
            else None
        )
        hole_area_frac = (
            float(
                (tissue_mask & (dist_hole > float(hole_distance_threshold))).sum()
                / n_tissue
            )
            if n_tissue > 0
            else None
        )

        if "is_blurred_gmm_2d" in df_grid_roi.columns:
            blurred_mask = df_grid_roi["is_blurred_gmm_2d"].to_numpy(bool, copy=False)
        elif "is_blurred_gmm" in df_grid_roi.columns:
            blurred_mask = df_grid_roi["is_blurred_gmm"].to_numpy(bool, copy=False)
        else:
            blurred_mask = np.zeros(len(df_grid_roi), dtype=bool)
        if "is_low_intensity" in df_grid_roi.columns:
            low_intensity_mask = df_grid_roi["is_low_intensity"].to_numpy(
                bool, copy=False
            )
        else:
            low_intensity_mask = np.zeros(len(df_grid_roi), dtype=bool)
        bad_mask = blurred_mask | low_intensity_mask

        usable_tissue_frac = (
            float((qc_tissue_mask & (~bad_mask)).sum() / n_qc_tissue)
            if n_qc_tissue > 0
            else None
        )

        # Largest contiguous bad zone (4-neighbour) over ROI grid.
        x_unique = np.unique(x1)
        y_unique = np.unique(y1)
        x_idx = np.searchsorted(x_unique, x1)
        y_idx = np.searchsorted(y_unique, y1)
        bad_grid = np.zeros((len(y_unique), len(x_unique)), dtype=bool)
        tissue_grid = np.zeros((len(y_unique), len(x_unique)), dtype=bool)
        bad_grid[y_idx, x_idx] = tissue_mask & bad_mask
        tissue_grid[y_idx, x_idx] = tissue_mask
        n_tissue_grid = int(tissue_grid.sum())
        cluster_zone_bad_fraction: float | None
        if n_tissue_grid > 0 and bad_grid.any():
            labels, n_comp = ndimage.label(
                bad_grid,
                structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8),
            )
            component_sizes = np.bincount(labels.ravel())
            largest_bad_component = int(component_sizes[1:].max()) if n_comp > 0 else 0
            cluster_zone_bad_fraction = float(largest_bad_component / n_tissue_grid)
        else:
            cluster_zone_bad_fraction = 0.0 if n_tissue_grid > 0 else None

        roi_metrics["morphology"] = {
            "edge_zone_frac": edge_zone_frac,
            "hole_area_frac": hole_area_frac,
            "usable_tissue_frac": usable_tissue_frac,
            "cluster_zone_bad_fraction": cluster_zone_bad_fraction,
            "edge_distance_threshold_px_ds": float(edge_distance_threshold),
            "hole_distance_threshold_px_ds": float(hole_distance_threshold),
            "min_tissue_coverage_for_qc": float(min_tissue_coverage_for_qc),
            "n_tissue_rois": n_tissue,
            "n_qc_tissue_rois": n_qc_tissue,
        }

    # Optional: Laplacian sharpness summary for section-6 report metrics.
    _qc = qc_thresholds or {}
    _ch_dapi = (_qc.get("channels") or {}).get("DAPI") or {}
    _focus_cut = _qc.get("focus") or {}
    lap_abs_floor_warn = _focus_cut.get("lap_var_absolute_floor_warn")
    lap_abs_floor_fail = _focus_cut.get("lap_var_absolute_floor_fail")

    if {"dapi_lap_var", "tissue_coverage"}.issubset(df_grid_roi.columns):
        lap_var = df_grid_roi["dapi_lap_var"].to_numpy(np.float64, copy=False)
        tissue_cov = df_grid_roi["tissue_coverage"].to_numpy(np.float64, copy=False)
        focus_col = (
            "dapi_focus_score"
            if "dapi_focus_score" in df_grid_roi.columns
            else ("focus_score" if "focus_score" in df_grid_roi.columns else None)
        )
        focus_vals = (
            df_grid_roi[focus_col].to_numpy(np.float64, copy=False)
            if focus_col is not None
            else np.full_like(lap_var, np.nan, dtype=np.float64)
        )

        mask = (
            (tissue_cov >= float(min_tissue_coverage_for_qc))
            & np.isfinite(lap_var)
            & np.isfinite(focus_vals)
            & (lap_var >= 0.0)
        )
        # P5: Exclude boundary ROIs (centre within window_size//2 of image edge)
        if "is_boundary_roi" in df_grid_roi.columns:
            boundary_mask = df_grid_roi["is_boundary_roi"].to_numpy(bool, copy=False)
            mask = mask & ~boundary_mask
            n_boundary_excluded = int(boundary_mask.sum())
        else:
            n_boundary_excluded = 0

        n_used = int(mask.sum())
        if n_used > 1:
            lap_used = lap_var[mask]
            focus_used = focus_vals[mask]

            med = float(np.nanmedian(lap_used))
            q25, q75 = np.nanpercentile(lap_used, [25, 75])
            iqr = float(q75 - q25)

            iqr_uniform = iqr <= 1e-12

            # P0: Absolute Laplacian variance floor (catches uniform defocus)
            lap_abs_status = "N/A"
            if isinstance(lap_abs_floor_fail, (int, float)) and np.isfinite(
                float(lap_abs_floor_fail)
            ):
                if med < float(lap_abs_floor_fail):
                    lap_abs_status = "FAIL"
                elif isinstance(lap_abs_floor_warn, (int, float)) and np.isfinite(
                    float(lap_abs_floor_warn)
                ):
                    if med < float(lap_abs_floor_warn):
                        lap_abs_status = "WARN"
                    else:
                        lap_abs_status = "PASS"
                else:
                    lap_abs_status = "PASS"
            elif isinstance(lap_abs_floor_warn, (int, float)) and np.isfinite(
                float(lap_abs_floor_warn)
            ):
                lap_abs_status = "WARN" if med < float(lap_abs_floor_warn) else "PASS"

            # P4: Guard Spearman correlation on uniform slides
            cv_focus = float(np.std(focus_used) / (np.mean(focus_used) + 1e-12))
            cv_lap = float(np.std(lap_used) / (np.mean(lap_used) + 1e-12))
            if cv_focus < 0.1 and cv_lap < 0.1:
                focus_lap_corr = None
                focus_lap_corr_status = "uniform"
            else:
                focus_lap_corr = float(
                    pd.Series(focus_used).corr(pd.Series(lap_used), method="spearman")
                )
                focus_lap_corr_status = None

            # Separation of 2D-GMM blur components in log1p(lap_var) space (Cohen's d).
            component_separation = None
            if "is_blurred_gmm_2d" in df_grid_roi.columns:
                blur2d = df_grid_roi["is_blurred_gmm_2d"].to_numpy(bool, copy=False)[
                    mask
                ]
                if blur2d.any() and (~blur2d).any():
                    lap_log = np.log1p(np.maximum(lap_used, 0.0))
                    lap_blur = lap_log[blur2d]
                    lap_focus = lap_log[~blur2d]
                    if lap_blur.size > 1 and lap_focus.size > 1:
                        n_b, n_f = lap_blur.size, lap_focus.size
                        v_blur = float(np.var(lap_blur, ddof=1))
                        v_focus = float(np.var(lap_focus, ddof=1))
                        pooled_sd = np.sqrt(
                            max(
                                ((n_b - 1) * v_blur + (n_f - 1) * v_focus)
                                / (n_b + n_f - 2),
                                0.0,
                            )
                            + 1e-12
                        )
                        component_separation = float(
                            abs(float(np.mean(lap_focus)) - float(np.mean(lap_blur)))
                            / pooled_sd
                        )

            roi_metrics["laplacian_sharpness"] = {
                "lap_sigma": float(lap_sigma) if lap_sigma is not None else None,
                "n_rois_used": n_used,
                "n_boundary_excluded": n_boundary_excluded,
                "min_tissue_coverage_for_qc": float(min_tissue_coverage_for_qc),
                "iqr_uniform": iqr_uniform,
                "lap_var_median_raw": med,
                "lap_var_absolute_floor_warn": float(lap_abs_floor_warn)
                if isinstance(lap_abs_floor_warn, (int, float))
                else None,
                "lap_var_absolute_floor_fail": float(lap_abs_floor_fail)
                if isinstance(lap_abs_floor_fail, (int, float))
                else None,
                "lap_var_absolute_status": lap_abs_status,
                "focus_lap_spearman_corr": focus_lap_corr,
                "focus_lap_corr_status": focus_lap_corr_status,
                "component_separation": component_separation,
            }

    # Save to JSON
    metrics_file = Path(outdir) / "roi_qc_metrics.json"
    with open(metrics_file, "w") as f:
        json.dump(roi_metrics, f, indent=2, default=str)
    logging.info(f"[OK] Saved tile QC metrics to {metrics_file}")


def save_pixel_focus_maps(focus_maps, outdir, compress=True):
    """
    Save per-pixel focus maps as compressed tiled TIFF files.

    Uses lz4 compression (5-10x faster than zlib) and parallel I/O via
    ThreadPoolExecutor to minimize wall-clock time on multi-map datasets.

    Args:
        focus_maps: dict returned by compute_all_focus_maps()
        outdir: Path to output directory
        compress: Use lz4 compression (default: True)
    """
    from concurrent.futures import ThreadPoolExecutor

    outdir = Path(outdir)
    items = [(name, arr) for name, arr in focus_maps.items() if arr is not None]

    if not items:
        logging.info("  No focus maps to save.")
        return

    def _save_one(name, array, outdir_path):
        t_start = time.time()
        path = outdir_path / f"{name}.tif"
        tile = (256, 256) if array.shape[0] >= 256 and array.shape[1] >= 256 else None
        tifffile.imwrite(
            str(path),
            array.astype(np.float32),
            bigtiff=True,
            compression="zstd" if compress else None,
            tile=tile,
        )
        elapsed = time.time() - t_start
        return name, path, array.shape, elapsed

    with ThreadPoolExecutor(max_workers=min(len(items), 4)) as executor:
        futures = [executor.submit(_save_one, name, arr, outdir) for name, arr in items]
        for future in futures:
            name, path, shape, elapsed = future.result()
            logging.info(f"  Saved {name}: {path} ({shape}, float32, {elapsed:.1f}s)")


def save_versions_file(outdir):
    """Save package versions to YAML file"""
    import matplotlib

    def get_version(package, package_name=None):
        """Safely get version of a package"""
        if package_name is None:
            package_name = (
                package.__name__ if hasattr(package, "__name__") else str(package)
            )

        try:
            if hasattr(package, "__version__"):
                return package.__version__
            else:
                # Try to get version via importlib
                import importlib.metadata

                return importlib.metadata.version(package_name)
        except Exception:
            return "unknown"

    # Get Python version
    python_version = (
        f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )

    # Get package versions safely
    packages = {
        "python": python_version,
        "numpy": get_version(np, "numpy"),
        "pandas": get_version(pd, "pandas"),
        "matplotlib": get_version(matplotlib, "matplotlib"),
        "seaborn": get_version(sns, "seaborn"),
        "click": get_version(click, "click"),
        "pathlib": "built-in",
    }

    # Get versions for other packages
    try:
        import skimage

        packages["scikit-image"] = get_version(skimage, "scikit-image")
    except Exception:
        packages["scikit-image"] = "unknown"

    try:
        packages["tifffile"] = get_version(tifffile, "tifffile")
    except Exception:
        packages["tifffile"] = "unknown"

    try:
        packages["zarr"] = get_version(zarr, "zarr")
    except Exception:
        packages["zarr"] = "unknown"

    try:
        import napari_skimage_regionprops

        packages["napari-skimage-regionprops"] = get_version(
            napari_skimage_regionprops, "napari-skimage-regionprops"
        )
    except Exception:
        packages["napari-skimage-regionprops"] = "unknown"

    try:
        packages["napari-simpleitk-image-processing"] = get_version(
            nsitk, "napari-simpleitk-image-processing"
        )
    except Exception:
        packages["napari-simpleitk-image-processing"] = "unknown"

    # Create YAML content
    yaml_content = f"""XENIUM_IMAGE_QC:
  python: {packages["python"]}
  numpy: {packages["numpy"]}
  pandas: {packages["pandas"]}
  matplotlib: {packages["matplotlib"]}
  seaborn: {packages["seaborn"]}
  scikit-image: {packages["scikit-image"]}
  tifffile: {packages["tifffile"]}
  zarr: {packages["zarr"]}
  napari-skimage-regionprops: {packages["napari-skimage-regionprops"]}
  napari-simpleitk-image-processing: {packages["napari-simpleitk-image-processing"]}
  click: {packages["click"]}
  pathlib: {packages["pathlib"]}
"""

    # Save to file
    versions_file = outdir / "versions.yml"
    with open(versions_file, "w") as f:
        f.write(yaml_content)

    logging.info(f"[OK] Saved package versions to {versions_file}")


# ===== FUNCTIONS FROM image_qc_mapping_to_cells.py =====


def load_spatial_data(data):
    """Load and prepare spatial data. This only works for the original Xenium bundle generated from XOA"""

    # Load cell masks
    cell_masks_zarr = open_zarr(data["cell_masks_path"])
    # cell_masks_zarr will only have one channel if its the result from xeniumranger import-segmentation
    cellseg_mask = np.array(cell_masks_zarr.get("masks").get("1"))

    # Load spatial data
    df_spatial = pd.read_parquet(
        data["cells_parquet_path"],
        columns=[
            "x_centroid",
            "y_centroid",
            "transcript_counts",
            "cell_area",
            "nucleus_area",
            "nucleus_count",
            "segmentation_method",
            "cell_id",
        ],
    )

    # Rename columns for simplicity
    df_spatial.rename(columns={"x_centroid": "x", "y_centroid": "y"}, inplace=True)

    # Convert to pixel coordinates
    df_spatial["x"] = (df_spatial.x / XENIUM_PIXEL_SIZE_UM).astype(int)
    df_spatial["y"] = (df_spatial.y / XENIUM_PIXEL_SIZE_UM).astype(int)

    # Measure CellID - OPTIMIZED: vectorized NumPy indexing (50-100x faster)
    x_coords = np.clip(df_spatial["x"].values, 0, cellseg_mask.shape[1] - 1)
    y_coords = np.clip(df_spatial["y"].values, 0, cellseg_mask.shape[0] - 1)
    df_spatial["CellID"] = cellseg_mask[y_coords, x_coords]

    # Load and merge clustering and UMAP data
    df_clusters = pd.read_csv(data["clusters_csv_path"]).rename(
        columns={"Barcode": "cell_id", "Cluster": "Cluster_kmeans10"}
    )
    df_UMAP = pd.read_csv(data["umap_path"]).rename(columns={"Barcode": "cell_id"})
    df_spatial = df_spatial.merge(df_clusters, on="cell_id").merge(
        df_UMAP, on="cell_id"
    )
    return df_spatial, cellseg_mask, cell_masks_zarr


def map_distances_to_cells(
    df_spatial,
    distance_map,
    distance_map2,
    dense_intensity_regions,
    downsample_factor=8,
):
    """
    Map distance maps and dense intensity region masks to cells based on cell centroid locations.

    Extracted from generate_sample_masks_and_distances() for cell-independent workflow.

    Parameters:
    -----------
    df_spatial : pandas DataFrame
        Cell DataFrame with 'x' and 'y' centroid columns
    distance_map : numpy.ndarray
        Distance to edge map (downsampled)
    distance_map2 : numpy.ndarray
        Distance to nearest hole map (downsampled)
    dense_intensity_regions : numpy.ndarray
        Dense intensity region mask (downsampled)
    downsample_factor : int, optional
        Downsampling factor (default: 8)

    Returns:
    --------
    pandas DataFrame
        df_spatial with added columns:
        - Distance-to-edge: Distance to sample edge for each cell
        - Distance-to-nearest-hole: Distance to nearest hole for each cell
        - In-Area-with-Dense-Intensity-Region: Binary indicator if cell is in dense intensity region area
    """
    # Vectorized coordinate mapping -- compute pixel indices once, clipped to array bounds
    X = np.clip(
        (df_spatial["x"].values / downsample_factor).astype(int),
        0,
        distance_map.shape[1] - 1,
    )
    Y = np.clip(
        (df_spatial["y"].values / downsample_factor).astype(int),
        0,
        distance_map.shape[0] - 1,
    )

    # Vectorized lookups -- three lines instead of three loops
    df_spatial["Distance-to-edge"] = distance_map[Y, X]
    df_spatial["Distance-to-nearest-hole"] = distance_map2[Y, X]
    df_spatial["Dense-Intensity-Region-ID"] = dense_intensity_regions[Y, X].astype(int)

    return df_spatial


def map_grid_roi_to_cells(df_grid_roi, df_cells, overlapping=False):
    """
    Map cell-independent grid tile focus scores to cells based on cell centroid location.

    This function replaces the previous cell-based tile approach with a cell-independent
    grid-based tile approach that provides uniform coverage across the tissue.

    For each cell, finds which tile its centroid (x, y) falls into and assigns
    that tile's focus scores to the cell.

    Parameters:
    -----------
    df_grid_roi : pandas DataFrame
        Grid ROI DataFrame with columns: roi_id, x1, x2, y1, y2, focus_score,
        focus_score_norm, raw_intensity, tissue_coverage, and channel-specific columns
    df_cells : pandas DataFrame
        Cell DataFrame with columns: x, y (centroid coordinates)
    overlapping : bool, optional
        Whether grid is overlapping. If True and cell falls into multiple ROIs,
        uses ROI with highest tissue_coverage (default: False)

    Returns:
    --------
    pandas DataFrame
        Cell DataFrame with added columns:
        - roi_id: Integer ID of ROI containing cell centroid (0, 1, 2, ...) or None if cell outside all ROIs
        - DAPI_mean_roi: Mean intensity from grid ROI (NaN if cell outside all ROIs)
        - DAPI_RFS_roi: Raw focus score from grid ROI (NaN if cell outside all ROIs)
        - DAPI_RFSnorm_roi: Normalized focus score from grid ROI (NaN if cell outside all ROIs)
        - roi_tissue_coverage: Tissue coverage of assigned ROI (NaN if cell outside all ROIs)
        - Additional channel-specific columns if available (boundary_*, intrna_*)
    """
    # Create output DataFrame
    df_cells_mapped = df_cells.copy()

    # Initialize columns using old naming convention for compatibility
    df_cells_mapped["roi_id"] = None
    df_cells_mapped["DAPI_mean_roi"] = np.nan
    df_cells_mapped["DAPI_RFS_roi"] = np.nan
    df_cells_mapped["DAPI_RFSnorm_roi"] = np.nan
    df_cells_mapped["roi_tissue_coverage"] = np.nan

    # Check for additional channel columns
    has_boundary = "boundary_focus_score" in df_grid_roi.columns
    has_intrna = "intrna_focus_score" in df_grid_roi.columns

    if has_boundary:
        df_cells_mapped["boundary_focus_score_roi"] = np.nan
        df_cells_mapped["boundary_focus_score_norm_roi"] = np.nan
        df_cells_mapped["boundary_intensity_roi"] = np.nan

    if has_intrna:
        df_cells_mapped["intrna_focus_score_roi"] = np.nan
        df_cells_mapped["intrna_focus_score_norm_roi"] = np.nan
        df_cells_mapped["intrna_intensity_roi"] = np.nan

    # OPTIMIZED: Fast vectorized lookup using 2D grid
    logging.info(f"  Mapping {len(df_cells):,} cells to {len(df_grid_roi):,} tiles...")

    # Get cell coordinates as arrays (faster than DataFrame access)
    cell_x = df_cells["x"].values
    cell_y = df_cells["y"].values
    n_cells = len(df_cells)

    # Create arrays to store matches
    cell_roi_matches = np.full(n_cells, -1, dtype=np.int32)  # -1 means no match
    cell_roi_tissue_coverage_arr = np.full(n_cells, np.nan, dtype=np.float64)

    # Always try the fast vectorized approach: build a 2D lookup grid
    is_regular_grid = False
    if len(df_grid_roi) > 1 and not overlapping:
        roi_size = df_grid_roi.iloc[0]["x2"] - df_grid_roi.iloc[0]["x1"]
        stride = roi_size  # non-overlapping grid

        x_min = df_grid_roi["x1"].min()
        y_min = df_grid_roi["y1"].min()
        x_max = df_grid_roi["x1"].max()
        y_max = df_grid_roi["y1"].max()

        n_cols = int(round((x_max - x_min) / stride)) + 1
        n_rows = int(round((y_max - y_min) / stride)) + 1

        # Build 2D grid: grid[row, col] -> index in df_grid_roi
        grid_lookup = np.full((n_rows, n_cols), -1, dtype=np.int32)
        roi_x1_arr = df_grid_roi["x1"].values
        roi_y1_arr = df_grid_roi["y1"].values
        gx_arr = np.round((roi_x1_arr - x_min) / stride).astype(int)
        gy_arr = np.round((roi_y1_arr - y_min) / stride).astype(int)

        # Check if this is a valid grid (no out-of-bounds)
        valid = (gx_arr >= 0) & (gx_arr < n_cols) & (gy_arr >= 0) & (gy_arr < n_rows)
        if valid.all():
            # Populate grid (vectorized)
            grid_lookup[gy_arr, gx_arr] = np.arange(len(df_grid_roi))

            # Vectorized lookup for all cells -- O(n_cells)
            cell_gx = np.clip(
                np.round((cell_x - x_min) / stride).astype(int), 0, n_cols - 1
            )
            cell_gy = np.clip(
                np.round((cell_y - y_min) / stride).astype(int), 0, n_rows - 1
            )
            cell_roi_idx = grid_lookup[cell_gy, cell_gx]  # vectorized!

            matched_mask = cell_roi_idx >= 0
            cell_roi_matches[matched_mask] = df_grid_roi["roi_id"].values[
                cell_roi_idx[matched_mask]
            ]
            cell_roi_tissue_coverage_arr[matched_mask] = df_grid_roi[
                "tissue_coverage"
            ].values[cell_roi_idx[matched_mask]]

            is_regular_grid = True
            logging.info(
                f"  Using fast grid lookup (stride={stride}px, grid={n_cols}x{n_rows})..."
            )

    if not is_regular_grid:
        # SLOW PATH: Irregular/filtered grid - fallback for grids that don't validate
        if len(df_grid_roi) > 1000:
            logging.info(
                f"  Processing {len(df_cells):,} cells against {len(df_grid_roi):,} tiles (irregular grid)..."
            )

        # Convert ROI boundaries to arrays for fast lookup
        roi_x1 = df_grid_roi["x1"].values
        roi_x2 = df_grid_roi["x2"].values
        roi_y1 = df_grid_roi["y1"].values
        roi_y2 = df_grid_roi["y2"].values
        roi_ids = df_grid_roi["roi_id"].values
        roi_tissue_coverage = df_grid_roi["tissue_coverage"].values

        for cell_idx in range(n_cells):
            x_cell = cell_x[cell_idx]
            y_cell = cell_y[cell_idx]

            matches = (
                (roi_x1 <= x_cell)
                & (x_cell < roi_x2)
                & (roi_y1 <= y_cell)
                & (y_cell < roi_y2)
            )

            if np.any(matches):
                match_idx = np.where(matches)[0][0]
                cell_roi_matches[cell_idx] = roi_ids[match_idx]
                cell_roi_tissue_coverage_arr[cell_idx] = roi_tissue_coverage[match_idx]

    # Check for GMM columns
    has_gmm_1d = "is_blurred_gmm" in df_grid_roi.columns
    has_gmm_2d = "is_blurred_gmm_2d" in df_grid_roi.columns
    has_blur_prob_1d = "blur_prob_gmm" in df_grid_roi.columns
    has_blur_prob_2d = "blur_prob_gmm_2d" in df_grid_roi.columns

    # Initialize GMM columns if available
    if has_gmm_1d:
        df_cells_mapped["is_blurred_gmm_roi"] = False
    if has_blur_prob_1d:
        df_cells_mapped["blur_prob_gmm_roi"] = np.nan
    if has_gmm_2d:
        df_cells_mapped["is_blurred_gmm_2d_roi"] = False
    if has_blur_prob_2d:
        df_cells_mapped["blur_prob_gmm_2d_roi"] = np.nan

    # Build vectorized column arrays from df_grid_roi for direct indexing
    # Map roi_id -> index in df_grid_roi for O(1) lookup
    roi_id_to_idx = pd.Series(
        range(len(df_grid_roi)), index=df_grid_roi["roi_id"].values
    )

    # Resolve column names (handle legacy naming)
    _dapi_intensity_col = (
        "dapi_intensity" if "dapi_intensity" in df_grid_roi.columns else "raw_intensity"
    )
    _dapi_focus_col = (
        "dapi_focus_score"
        if "dapi_focus_score" in df_grid_roi.columns
        else "focus_score"
    )
    _dapi_focus_norm_col = (
        "dapi_focus_score_norm"
        if "dapi_focus_score_norm" in df_grid_roi.columns
        else "focus_score_norm"
    )

    dapi_intensity_arr = df_grid_roi[_dapi_intensity_col].values
    dapi_focus_arr = df_grid_roi[_dapi_focus_col].values
    dapi_focus_norm_arr = df_grid_roi[_dapi_focus_norm_col].values

    # Assign ROI values to cells using vectorized operations
    matched_mask = cell_roi_matches >= 0
    matched_roi_ids = cell_roi_matches[matched_mask]

    # Map matched roi_ids to df_grid_roi indices for vectorized column access
    matched_df_indices = roi_id_to_idx[matched_roi_ids].values

    # Assign ROI ID and tissue coverage
    df_cells_mapped.loc[matched_mask, "roi_id"] = matched_roi_ids
    df_cells_mapped.loc[matched_mask, "roi_tissue_coverage"] = (
        cell_roi_tissue_coverage_arr[matched_mask]
    )

    # Assign DAPI data via direct array indexing (no dict lookups)
    df_cells_mapped.loc[matched_mask, "DAPI_mean_roi"] = dapi_intensity_arr[
        matched_df_indices
    ]
    df_cells_mapped.loc[matched_mask, "DAPI_RFS_roi"] = dapi_focus_arr[
        matched_df_indices
    ]
    df_cells_mapped.loc[matched_mask, "DAPI_RFSnorm_roi"] = dapi_focus_norm_arr[
        matched_df_indices
    ]

    # Assign GMM classifications if available
    if has_gmm_1d:
        df_cells_mapped.loc[matched_mask, "is_blurred_gmm_roi"] = df_grid_roi[
            "is_blurred_gmm"
        ].values[matched_df_indices]
    if has_blur_prob_1d:
        df_cells_mapped.loc[matched_mask, "blur_prob_gmm_roi"] = df_grid_roi[
            "blur_prob_gmm"
        ].values[matched_df_indices]
    if has_gmm_2d:
        df_cells_mapped.loc[matched_mask, "is_blurred_gmm_2d_roi"] = df_grid_roi[
            "is_blurred_gmm_2d"
        ].values[matched_df_indices]
    if has_blur_prob_2d:
        df_cells_mapped.loc[matched_mask, "blur_prob_gmm_2d_roi"] = df_grid_roi[
            "blur_prob_gmm_2d"
        ].values[matched_df_indices]

    # Assign additional channel data if available
    if has_boundary:
        df_cells_mapped.loc[matched_mask, "boundary_focus_score_roi"] = df_grid_roi[
            "boundary_focus_score"
        ].values[matched_df_indices]
        df_cells_mapped.loc[matched_mask, "boundary_focus_score_norm_roi"] = (
            df_grid_roi["boundary_focus_score_norm"].values[matched_df_indices]
        )
        df_cells_mapped.loc[matched_mask, "boundary_intensity_roi"] = df_grid_roi[
            "boundary_intensity"
        ].values[matched_df_indices]

    if has_intrna:
        df_cells_mapped.loc[matched_mask, "intrna_focus_score_roi"] = df_grid_roi[
            "intrna_focus_score"
        ].values[matched_df_indices]
        df_cells_mapped.loc[matched_mask, "intrna_focus_score_norm_roi"] = df_grid_roi[
            "intrna_focus_score_norm"
        ].values[matched_df_indices]
        df_cells_mapped.loc[matched_mask, "intrna_intensity_roi"] = df_grid_roi[
            "intrna_intensity"
        ].values[matched_df_indices]

    # Handle overlapping grids: if cell falls into multiple ROIs, use highest tissue_coverage
    if overlapping:
        # Find cells with multiple ROI assignments
        cell_roi_counts = df_cells_mapped.groupby(df_cells_mapped.index)[
            "roi_id"
        ].count()
        cells_with_multiple = cell_roi_counts[cell_roi_counts > 1].index

        if len(cells_with_multiple) > 0:
            # For each cell with multiple ROIs, find the one with highest tissue_coverage
            for cell_idx in cells_with_multiple:
                # Find all ROIs this cell falls into
                x_cell = df_cells.loc[cell_idx, "x"]
                y_cell = df_cells.loc[cell_idx, "y"]

                matching_rois = df_grid_roi[
                    (df_grid_roi["x1"] <= x_cell)
                    & (x_cell < df_grid_roi["x2"])
                    & (df_grid_roi["y1"] <= y_cell)
                    & (y_cell < df_grid_roi["y2"])
                ]

                if len(matching_rois) > 0:
                    # Select ROI with highest tissue_coverage
                    best_roi = matching_rois.loc[
                        matching_rois["tissue_coverage"].idxmax()
                    ]

                    # Update cell with best ROI
                    df_cells_mapped.loc[cell_idx, "roi_id"] = best_roi["roi_id"]
                    df_cells_mapped.loc[cell_idx, "DAPI_mean_roi"] = best_roi.get(
                        "dapi_intensity", best_roi.get("raw_intensity", np.nan)
                    )
                    df_cells_mapped.loc[cell_idx, "DAPI_RFS_roi"] = best_roi.get(
                        "dapi_focus_score", best_roi.get("focus_score", np.nan)
                    )
                    df_cells_mapped.loc[cell_idx, "DAPI_RFSnorm_roi"] = best_roi.get(
                        "dapi_focus_score_norm",
                        best_roi.get("focus_score_norm", np.nan),
                    )
                    df_cells_mapped.loc[cell_idx, "roi_tissue_coverage"] = best_roi[
                        "tissue_coverage"
                    ]

                    # Update GMM classifications if available
                    if has_gmm_1d:
                        df_cells_mapped.loc[cell_idx, "is_blurred_gmm_roi"] = (
                            best_roi.get("is_blurred_gmm", False)
                        )
                    if has_blur_prob_1d:
                        df_cells_mapped.loc[cell_idx, "blur_prob_gmm_roi"] = (
                            best_roi.get("blur_prob_gmm", np.nan)
                        )
                    if has_gmm_2d:
                        df_cells_mapped.loc[cell_idx, "is_blurred_gmm_2d_roi"] = (
                            best_roi.get("is_blurred_gmm_2d", False)
                        )
                    if has_blur_prob_2d:
                        df_cells_mapped.loc[cell_idx, "blur_prob_gmm_2d_roi"] = (
                            best_roi.get("blur_prob_gmm_2d", np.nan)
                        )

                    if has_boundary:
                        df_cells_mapped.loc[cell_idx, "boundary_focus_score_roi"] = (
                            best_roi.get("boundary_focus_score", np.nan)
                        )
                        df_cells_mapped.loc[
                            cell_idx, "boundary_focus_score_norm_roi"
                        ] = best_roi.get("boundary_focus_score_norm", np.nan)
                        df_cells_mapped.loc[cell_idx, "boundary_intensity_roi"] = (
                            best_roi.get("boundary_intensity", np.nan)
                        )

                    if has_intrna:
                        df_cells_mapped.loc[cell_idx, "intrna_focus_score_roi"] = (
                            best_roi.get("intrna_focus_score", np.nan)
                        )
                        df_cells_mapped.loc[cell_idx, "intrna_focus_score_norm_roi"] = (
                            best_roi.get("intrna_focus_score_norm", np.nan)
                        )
                        df_cells_mapped.loc[cell_idx, "intrna_intensity_roi"] = (
                            best_roi.get("intrna_intensity", np.nan)
                        )

    return df_cells_mapped


def load_roi_blur_threshold(outdir):
    """
    Load ROI blur threshold configuration from JSON file saved by image_qc_roi_processing.py.

    Parameters:
    -----------
    outdir : Path
        Output directory where roi_blur_threshold.json should be located

    Returns:
    --------
    tuple
        (roi_focus_score_threshold, roi_intensity_threshold) or (None, None) if not found
    """
    threshold_json = Path(outdir) / "roi_blur_threshold.json"
    if not threshold_json.exists():
        return None, None

    try:
        with open(threshold_json, "r") as f:
            threshold_config = json.load(f)
        roi_threshold = threshold_config.get("roi_focus_score_threshold")
        intensity_threshold = threshold_config.get("roi_intensity_threshold")
        return roi_threshold, intensity_threshold
    except (json.JSONDecodeError, KeyError) as e:
        logging.warning(f"  Warning: Could not load threshold configuration: {e}")
        return None, None


def create_final_merged_data(
    df_spatial,
    myData,
    roi_data=None,
    ccfs_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
    edge_distance_threshold=-25,
    hole_distance_threshold=-25,
    roi_threshold=None,
    roi_intensity_threshold=None,
):
    """Create final merged dataset with boolean columns for thresholds"""

    # Link the Spatial Cell matrix (df_spatial) with the FocusScore table
    new_df = pd.merge(df_spatial, myData, how="inner", on="CellID")

    # Create new segmentation names with colour palette for plotting.
    # XR import-segmentation (e.g. from the CROP subworkflow) writes
    # `segmentation_method = 'Imported Cell Segmentation'` into cells.parquet,
    # which doesn't substring-match any of OBA's three categories
    # ("boundary"/"interior"/"nucleus"). Without this 4th row the
    # match-filter below drops every cell and downstream UMAP / cluster /
    # nuclear-texture / GMM figures fail with "no numeric data to plot"
    # or matplotlib Normalize errors on empty arrays.
    seg = pd.DataFrame(
        {
            "segmentation": [
                "boundary",
                "interior",
                "nucleus",
                "Imported Cell Segmentation",
            ],
            "segPal": ["#FFABC3", "#A9A800", "#A9CEFF", "#7FBFFF"],
            "join": [1, 1, 1, 1],
        }
    )

    # Link the Spatial Cell matrix to the new segmentation colour palette
    new_df["join"] = 1
    new_df = new_df.merge(seg, on="join").drop("join", axis=1)
    seg.drop("join", axis=1, inplace=True)
    new_df["match"] = new_df.apply(
        lambda x: x.segmentation_method.find(x.segmentation), axis=1
    ).ge(0)
    new_df = new_df[new_df["match"]]

    # Add boolean columns for various thresholds (nuclei-based method)
    new_df["is_low_nuclear_texture"] = new_df["CCFS_DAPI"] <= ccfs_threshold
    new_df["is_high_nuclear_texture"] = new_df["CCFS_DAPI"] > ccfs_threshold
    new_df["is_near_edge"] = new_df["Distance-to-edge"] > edge_distance_threshold
    new_df["is_far_from_edge"] = new_df["Distance-to-edge"] <= edge_distance_threshold
    new_df["is_near_hole"] = (
        new_df["Distance-to-nearest-hole"] > hole_distance_threshold
    )
    new_df["is_far_from_hole"] = (
        new_df["Distance-to-nearest-hole"] <= hole_distance_threshold
    )
    # Create boolean indicator for cells in any dense intensity region (backward compatibility)
    new_df["has_dense_intensity_regions"] = new_df["Dense-Intensity-Region-ID"] > 0

    # Merge ROI data if provided
    calculated_roi_threshold = None
    if roi_data is not None:
        # Merge ROI data using x and y coordinates
        new_df = pd.merge(
            new_df, roi_data, on=["x", "y"], how="left", suffixes=("", "_roi_merge")
        )

        # Use new threshold approach: raw score percentile + intensity threshold
        # If roi_threshold is provided (e.g., for backward compatibility), use it
        # Otherwise, calculate from raw scores using configured parameters
        if roi_threshold is None:
            # roi_threshold should have been calculated from df_grid_roi and passed here
            # If not provided, we can't calculate it here (need df_grid_roi)
            # This should not happen in normal flow, but provide fallback
            logging.warning(
                "  Warning: roi_threshold not provided, using default calculation"
            )
            roi_threshold = -1.0  # Old default for backward compatibility

        # Use intensity threshold from configuration if not provided
        # If None, it means threshold config wasn't loaded (backward compatibility)
        # In that case, skip intensity check and only use focus score threshold
        if roi_intensity_threshold is None:
            roi_intensity_threshold = (
                None  # Will skip intensity check in classification
            )

        calculated_roi_threshold = roi_threshold

        # Add boolean columns for tile-based blur detection
        # New approach: blurred if (raw_focus_score <= threshold) OR (intensity < intensity_threshold)
        # Use raw scores (DAPI_RFS_roi) not normalized (DAPI_RFSnorm_roi)
        has_intensity = "DAPI_mean_roi" in new_df.columns
        has_raw_focus = "DAPI_RFS_roi" in new_df.columns

        if has_raw_focus:
            # Combined threshold: focus score OR intensity (if intensity threshold provided)
            if has_intensity and roi_intensity_threshold is not None:
                new_df["is_blurred_roi"] = (new_df["DAPI_RFS_roi"] <= roi_threshold) | (
                    new_df["DAPI_mean_roi"] < roi_intensity_threshold
                )
            else:
                # Fallback: just use focus score if intensity not available or threshold not provided
                new_df["is_blurred_roi"] = new_df["DAPI_RFS_roi"] <= roi_threshold
        else:
            # Fallback: use normalized scores if raw scores not available (backward compatibility)
            new_df["is_blurred_roi"] = new_df["DAPI_RFSnorm_roi"] <= roi_threshold

        new_df["is_high_focus_roi"] = ~new_df["is_blurred_roi"]

    return new_df, calculated_roi_threshold


# ===== CELL-LEVEL PLOTTING FUNCTIONS =====


def plot_nuclear_texture_proportions(
    new_df,
    figures_dir,
    figures_source_dir,
    texture_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
    GROUP_BY_COLUMN="Cluster_kmeans10",
):
    """
    Create a stacked bar plot showing proportion of cells with high and low blur scores per cluster.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing CCFS_DAPI and Cluster_kmeans10 columns
    texture_threshold : float, optional
        Threshold to classify cells as high/low nuclear texture (default: DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD)
    """
    # Create nuclear texture categories (create temporary column to avoid modifying original)
    df_temp = new_df.copy()
    # Low CCFS_DAPI values (<= threshold) = low nuclear texture quality
    # High CCFS_DAPI values (> threshold) = high nuclear texture quality
    df_temp["texture_category"] = np.where(
        df_temp["CCFS_DAPI"] <= texture_threshold, "Low Quality", "High Quality"
    )

    # Calculate proportions
    proportions = (
        df_temp.groupby([GROUP_BY_COLUMN, "texture_category"])
        .size()
        .unstack(fill_value=0)
    )
    proportions = (
        proportions.div(proportions.sum(axis=1), axis=0) * 100
    )  # Convert to percentages

    # Set up the plot
    fig = plt.figure(figsize=(12, 6))

    # Create stacked bar plot with specified colors
    # Colors are applied in alphabetical order: 'High Quality' (blue), 'Low Quality' (red)
    proportions.plot(
        kind="bar",
        stacked=True,
        color=["#1f77b4", "#d62728"],  # Blue for high quality, red for low quality
        figsize=(12, 6),
    )

    # Customize the plot
    plt.title(
        f"Proportion of Cells by Nuclear Texture Quality (Threshold: {texture_threshold})",
        fontsize=14,
        pad=20,
    )
    plt.xlabel(GROUP_BY_COLUMN, fontsize=12)
    plt.ylabel("Percentage of Cells", fontsize=12)

    # Rotate x-axis labels if needed
    plt.xticks(rotation=45, ha="right")

    # Add percentage labels on the bars
    for c in plt.gca().containers:
        # Add labels
        plt.gca().bar_label(c, fmt="%.1f%%", label_type="center")

    # Add a grid for better readability
    plt.grid(True, axis="y", linestyle="--", alpha=0.7)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save figure instead of showing
    plt.savefig(
        figures_dir / "nuclear_texture_proportions.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV
    df_texture_proportions = (
        df_temp.groupby([GROUP_BY_COLUMN, "texture_category"])
        .size()
        .reset_index(name="count")
    )
    df_texture_proportions.to_csv(
        figures_source_dir / "nuclear_texture_proportions.csv", index=False
    )

    # Print summary statistics
    logging.info("Summary statistics by cluster:")
    summary_stats = (
        df_temp.groupby([GROUP_BY_COLUMN, "texture_category"])
        .size()
        .unstack(fill_value=0)
    )
    logging.info("Count of cells by texture category:")
    logging.info(summary_stats)
    logging.info("Percentage of cells by texture category:")
    logging.info(summary_stats.div(summary_stats.sum(axis=1), axis=0) * 100)


def plot_tile_blur_proportions_roi(
    new_df,
    figures_dir,
    figures_source_dir,
    roi_threshold=-1.0,
    GROUP_BY_COLUMN="Cluster_kmeans10",
):
    """
    Create a stacked bar plot showing proportion of cells with high and low tile-based blur scores per cluster.
    Uses GMM 2D classification if available, otherwise falls back to threshold-based method.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing DAPI_RFSnorm_roi, is_blurred_roi (or is_blurred_gmm_2d_roi), and Cluster_kmeans10 columns
    roi_threshold : float, optional
        Threshold to classify cells as high/low blur (default: -1.0) - used only if GMM 2D not available
    """
    # Check if ROI data is available - prefer GMM 2D, fall back to threshold-based
    has_gmm_2d = "is_blurred_gmm_2d_roi" in new_df.columns
    has_threshold = "is_blurred_roi" in new_df.columns

    if not has_gmm_2d and not has_threshold:
        logging.warning(
            "Warning: Tile-based blur scores not available. Skipping tile blur score proportions plot."
        )
        return

    # Create blur score categories (create temporary column to avoid modifying original)
    df_temp = new_df.copy()

    # Use GMM 2D if available, otherwise use threshold-based
    if has_gmm_2d:
        df_temp["blur_category"] = np.where(
            df_temp["is_blurred_gmm_2d_roi"], "Blurred", "In Focus"
        )
        method_name = "2D GMM"
    else:
        df_temp["blur_category"] = np.where(
            df_temp["is_blurred_roi"], "Blurred", "In Focus"
        )
        method_name = "Threshold"

    # Calculate proportions
    proportions = (
        df_temp.groupby([GROUP_BY_COLUMN, "blur_category"]).size().unstack(fill_value=0)
    )
    proportions = (
        proportions.div(proportions.sum(axis=1), axis=0) * 100
    )  # Convert to percentages

    # Set up the plot
    fig = plt.figure(figsize=(12, 6))

    # Create stacked bar plot with specified colors
    # Colors are applied in alphabetical order: 'Blurred' (red), 'In Focus' (blue)
    proportions.plot(
        kind="bar",
        stacked=True,
        color=["#d62728", "#1f77b4"],  # Red for blurred, blue for in focus
        figsize=(12, 6),
    )

    # Customize the plot
    if has_gmm_2d:
        plt.title(
            "Proportion of Cells by Tile-based Blur Score (2D GMM Classification)",
            fontsize=14,
            pad=20,
        )
    else:
        # Note: roi_threshold is in raw score units, classification uses combined approach
        # Try to get intensity threshold from the threshold config if available
        intensity_threshold_display = getattr(
            plot_tile_blur_proportions_roi, "_intensity_threshold", None
        )
        if intensity_threshold_display is None:
            intensity_threshold_display = 20.0  # Default
        plt.title(
            f"Proportion of Cells by Tile-based Blur Score\n(Raw score threshold: {roi_threshold:.2f}, Intensity threshold: {intensity_threshold_display})",
            fontsize=14,
            pad=20,
        )
    plt.xlabel(GROUP_BY_COLUMN, fontsize=12)
    plt.ylabel("Percentage of Cells", fontsize=12)

    # Rotate x-axis labels if needed
    plt.xticks(rotation=45, ha="right")

    # Add percentage labels on the bars
    for c in plt.gca().containers:
        # Add labels
        plt.gca().bar_label(c, fmt="%.1f%%", label_type="center")

    # Add a grid for better readability
    plt.grid(True, axis="y", linestyle="--", alpha=0.7)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save figure instead of showing
    plt.savefig(
        figures_dir / "tile_blur_proportions_roi.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV
    df_blur_proportions = (
        df_temp.groupby([GROUP_BY_COLUMN, "blur_category"])
        .size()
        .reset_index(name="count")
    )
    df_blur_proportions.to_csv(
        figures_source_dir / "tile_blur_proportions_roi.csv", index=False
    )

    # Print summary statistics
    logging.info(f"Summary statistics by cluster (tile-based, {method_name}):")
    summary_stats = (
        df_temp.groupby([GROUP_BY_COLUMN, "blur_category"]).size().unstack(fill_value=0)
    )
    logging.info("Count of cells by blur category:")
    logging.info(summary_stats)
    logging.info("Percentage of cells by blur category:")
    logging.info(summary_stats.div(summary_stats.sum(axis=1), axis=0) * 100)


def plot_cell_focus_distribution(
    new_df,
    figures_dir,
    figures_source_dir,
    GROUP_BY_COLUMN="Cluster_kmeans10",
):
    """Plot cell-level focus score distributions using tile-propagated metrics.

    Creates a 2×2 panel:
      - Top-left: histogram of DAPI_RFSnorm_roi (tile focus score per cell)
      - Top-right: histogram coloured by GMM blur classification
      - Bottom-left: density of DAPI_RFSnorm_roi per cluster
      - Bottom-right: density of blur_prob_gmm_2d_roi per cluster
    """
    has_focus = "DAPI_RFSnorm_roi" in new_df.columns
    has_gmm_prob = "blur_prob_gmm_2d_roi" in new_df.columns
    has_gmm_class = "is_blurred_gmm_2d_roi" in new_df.columns

    if not has_focus:
        logging.warning(
            "No DAPI_RFSnorm_roi column — skipping cell focus distribution plot"
        )
        return

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    focus_vals = new_df["DAPI_RFSnorm_roi"].dropna()

    # --- Top-left: overall focus score histogram ---
    ax = axes[0, 0]
    ax.hist(focus_vals, bins=60, alpha=0.7, edgecolor="black", color="steelblue")
    ax.set_xlabel("Tile focus score (DAPI_RFSnorm_roi)", fontsize=11)
    ax.set_ylabel("Number of cells", fontsize=11)
    ax.set_title("Cell-level focus score distribution", fontsize=13)
    stats_text = f"n={len(focus_vals):,}\nMedian={focus_vals.median():.4f}\nMean={focus_vals.mean():.4f}"
    ax.text(
        0.95,
        0.95,
        stats_text,
        transform=ax.transAxes,
        va="top",
        ha="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
        fontsize=10,
    )
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # --- Top-right: histogram split by GMM blur class ---
    ax = axes[0, 1]
    if has_gmm_class:
        gmm_col = new_df["is_blurred_gmm_2d_roi"].astype(bool)
        valid = new_df["DAPI_RFSnorm_roi"].notna()
        blurred = new_df.loc[valid & gmm_col, "DAPI_RFSnorm_roi"]
        sharp = new_df.loc[valid & ~gmm_col, "DAPI_RFSnorm_roi"]
        ax.hist(
            sharp,
            bins=60,
            alpha=0.6,
            label=f"In focus ({len(sharp):,})",
            color="#1f77b4",
        )
        ax.hist(
            blurred,
            bins=60,
            alpha=0.6,
            label=f"Blurred ({len(blurred):,})",
            color="#d62728",
        )
        ax.legend(fontsize=10)
        ax.set_title("Focus score by GMM blur classification", fontsize=13)
    else:
        ax.hist(focus_vals, bins=60, alpha=0.7, color="steelblue")
        ax.set_title("Focus score distribution (no GMM data)", fontsize=13)
    ax.set_xlabel("Tile focus score (DAPI_RFSnorm_roi)", fontsize=11)
    ax.set_ylabel("Number of cells", fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # --- Bottom-left: focus score density per cluster ---
    ax = axes[1, 0]
    if GROUP_BY_COLUMN in new_df.columns:
        clusters = sorted(new_df[GROUP_BY_COLUMN].dropna().unique())
        for cl in clusters:
            vals = new_df.loc[
                new_df[GROUP_BY_COLUMN] == cl, "DAPI_RFSnorm_roi"
            ].dropna()
            if len(vals) > 10:
                vals.plot.kde(ax=ax, label=f"C{cl}", alpha=0.7)
        ax.legend(fontsize=8, ncol=2)
    else:
        focus_vals.plot.kde(ax=ax, color="steelblue")
    ax.set_xlabel("Tile focus score (DAPI_RFSnorm_roi)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title("Focus score density by cluster", fontsize=13)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    # --- Bottom-right: GMM blur probability density per cluster ---
    ax = axes[1, 1]
    if has_gmm_prob and GROUP_BY_COLUMN in new_df.columns:
        clusters = sorted(new_df[GROUP_BY_COLUMN].dropna().unique())
        for cl in clusters:
            vals = new_df.loc[
                new_df[GROUP_BY_COLUMN] == cl, "blur_prob_gmm_2d_roi"
            ].dropna()
            if len(vals) > 10:
                vals.plot.kde(ax=ax, label=f"C{cl}", alpha=0.7)
        ax.axvline(
            x=0.5, color="red", linestyle="--", alpha=0.6, label="Blur threshold (0.5)"
        )
        ax.legend(fontsize=8, ncol=2)
        ax.set_title("GMM blur probability by cluster", fontsize=13)
    elif has_gmm_prob:
        new_df["blur_prob_gmm_2d_roi"].dropna().plot.kde(ax=ax, color="steelblue")
        ax.axvline(x=0.5, color="red", linestyle="--", alpha=0.6)
        ax.set_title("GMM blur probability distribution", fontsize=13)
    else:
        ax.text(
            0.5,
            0.5,
            "No GMM blur probability data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="gray",
        )
        ax.set_title("GMM blur probability (not available)", fontsize=13)
    ax.set_xlabel("P(blur component)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.savefig(
        figures_dir / "cell_focus_distribution.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save source data
    src_cols = ["DAPI_RFSnorm_roi"]
    if has_gmm_class:
        src_cols.append("is_blurred_gmm_2d_roi")
    if has_gmm_prob:
        src_cols.append("blur_prob_gmm_2d_roi")
    if GROUP_BY_COLUMN in new_df.columns:
        src_cols.append(GROUP_BY_COLUMN)
    src = new_df[[c for c in src_cols if c in new_df.columns]].dropna(
        subset=["DAPI_RFSnorm_roi"]
    )
    src.to_csv(figures_source_dir / "cell_focus_distribution.csv", index=False)


def plot_nuclear_texture_density(
    new_df,
    figures_dir,
    figures_source_dir,
    GROUP_BY_COLUMN="Cluster_kmeans10",
    ccfs_low_texture_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
):
    """
    Create a density plot of CCFS_DAPI values grouped by Cluster_kmeans10.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing CCFS_DAPI and Cluster_kmeans10 columns
    """
    # Set up the plot
    fig = plt.figure(figsize=(12, 6))

    # Create the density plot
    ax = sns.kdeplot(
        data=new_df,
        x="CCFS_DAPI",
        hue=GROUP_BY_COLUMN,
        palette="husl",
        common_norm=False,  # Normalize each cluster separately
        fill=True,
        alpha=0.3,
    )  # Make the fill semi-transparent

    # Add a vertical line for the blur threshold
    plt.axvline(
        x=ccfs_low_texture_threshold,
        color="black",
        linestyle="--",
        alpha=0.5,
        label=f"Low Texture Threshold ({ccfs_low_texture_threshold})",
    )

    # Customize the plot
    plt.title(
        f"Distribution of Nuclear Texture Scores (CCFS_DAPI) by {GROUP_BY_COLUMN}",
        fontsize=14,
        pad=20,
    )
    plt.xlabel("CCFS_DAPI (Nuclear Texture Score)", fontsize=12)
    plt.ylabel("Density", fontsize=12)

    # Add a grid for better readability
    plt.grid(True, linestyle="--", alpha=0.7)

    # Get the current legend
    legend = ax.get_legend()

    # Get the handles and labels
    handles = legend.legend_handles
    labels = [f"Cluster {label.get_text()}" for label in legend.get_texts()]

    # Add the threshold line to the legend
    threshold_line = plt.Line2D([0], [0], color="black", linestyle="--", alpha=0.5)
    handles = [threshold_line] + handles
    labels = ["Default Threshold"] + labels

    # Create new legend
    plt.legend(
        handles,
        labels,
        title=GROUP_BY_COLUMN,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save figure instead of showing
    plt.savefig(
        figures_dir / "nuclear_texture_density.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV
    df_density_data = new_df[["CCFS_DAPI", GROUP_BY_COLUMN]].copy()
    df_density_data.to_csv(
        figures_source_dir / "nuclear_texture_density.csv", index=False
    )


def plot_nuclear_texture_vs_transcripts(
    new_df,
    figures_dir,
    figures_source_dir,
    log_scale=False,
    ccfs_low_texture_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
):
    """
    Create a scatter plot of nuclear texture scores (CCFS_DAPI) vs transcript counts.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing CCFS_DAPI and transcript_counts columns
    log_scale : bool, optional
        Whether to use log scale for x-axis (default: False)
    """
    # Set up the plot
    fig = plt.figure(figsize=(12, 6))

    # Create scatter plot
    plt.scatter(
        new_df["transcript_counts"],
        new_df["CCFS_DAPI"],
        alpha=0.5,  # Make points semi-transparent
        s=10,  # Set point size
        rasterized=True,
    )

    # Add a horizontal line for the blur threshold
    plt.axhline(
        y=ccfs_low_texture_threshold,
        color="red",
        linestyle="--",
        alpha=0.5,
        label=f"Low Texture Threshold ({ccfs_low_texture_threshold})",
    )

    # Customize the plot
    title = "Relationship between Nuclear Texture Score and Transcript Counts"
    if log_scale:
        title += " (Log Scale)"
    plt.title(title, fontsize=14, pad=20)
    plt.ylabel("CCFS_DAPI (Nuclear Texture Score)", fontsize=12)
    plt.xlabel("Transcript Counts", fontsize=12)

    # Set log scale for both axes if requested
    if log_scale:
        plt.xscale("log")
        plt.xlabel("Transcript Counts (log scale)", fontsize=12)
        plt.yscale("log")
        plt.ylabel("CCFS_DAPI (Nuclear Texture Score) (log scale)", fontsize=12)

    # Add a grid for better readability
    plt.grid(True, linestyle="--", alpha=0.7)

    # Add legend
    plt.legend()

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save figure instead of showing
    suffix = "_log" if log_scale else ""
    plt.savefig(
        figures_dir / f"nuclear_texture_vs_transcripts{suffix}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data as CSV
    df_scatter_data = new_df[["transcript_counts", "CCFS_DAPI"]].copy()
    df_scatter_data.to_csv(
        figures_source_dir / f"nuclear_texture_vs_transcripts{suffix}.csv", index=False
    )

    # Print summary statistics
    logging.info("Summary statistics:")
    logging.info("Nuclear Texture Score (CCFS_DAPI):")
    logging.info(new_df["CCFS_DAPI"].describe().round(4))
    logging.info("Transcript Counts:")
    logging.info(new_df["transcript_counts"].describe())


def plot_gmm_focus_vs_transcripts(
    new_df,
    figures_dir,
    figures_source_dir,
):
    """Scatter plot of tile-level GMM focus score vs transcript counts, coloured by blur class.

    Mirrors plot_nuclear_texture_vs_transcripts but uses the tile-level Laplacian
    GMM focus score (DAPI_RFSnorm_roi) propagated to cells, with points coloured
    by GMM blur classification.
    """
    focus_col = "DAPI_RFSnorm_roi"
    tx_col = "transcript_counts"
    gmm_col = "is_blurred_gmm_2d_roi"

    if focus_col not in new_df.columns or tx_col not in new_df.columns:
        logging.warning(
            "Missing %s or %s — skipping GMM focus vs transcripts plot.",
            focus_col,
            tx_col,
        )
        return

    df = new_df[
        [tx_col, focus_col] + ([gmm_col] if gmm_col in new_df.columns else [])
    ].dropna()
    if df.empty:
        return

    has_gmm = gmm_col in df.columns
    fig, ax = plt.subplots(figsize=(12, 6))

    if has_gmm:
        sharp = df[~df[gmm_col].astype(bool)]
        blurred = df[df[gmm_col].astype(bool)]
        ax.scatter(
            sharp[tx_col],
            sharp[focus_col],
            alpha=0.3,
            s=8,
            color="#1f77b4",
            label=f"In Focus ({len(sharp):,})",
            rasterized=True,
        )
        ax.scatter(
            blurred[tx_col],
            blurred[focus_col],
            alpha=0.3,
            s=8,
            color="#d62728",
            label=f"Blurred ({len(blurred):,})",
            rasterized=True,
        )
        ax.legend(fontsize=10)
    else:
        ax.scatter(df[tx_col], df[focus_col], alpha=0.3, s=8, rasterized=True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Transcript Counts (log scale)", fontsize=12)
    ax.set_ylabel("Tile Focus Score (DAPI_RFSnorm_roi, log scale)", fontsize=12)
    ax.set_title("GMM Tile Focus Score vs Transcript Counts", fontsize=14, pad=20)
    ax.grid(True, linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.savefig(
        figures_dir / "gmm_focus_vs_transcripts.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save source data
    src_cols = [tx_col, focus_col]
    if has_gmm:
        src_cols.append(gmm_col)
    df[src_cols].to_csv(
        figures_source_dir / "gmm_focus_vs_transcripts.csv", index=False
    )

    logging.info("GMM focus vs transcripts: %d cells plotted.", len(df))


def plot_ccfs_vs_roi_comparison(new_df, figures_dir, figures_source_dir):
    """
    Create comparison plots between CCFS (nuclei-based) and cell-independent tile-based focus scores.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing both CCFS_DAPI and DAPI_RFSnorm_roi columns
        (DAPI_RFSnorm_roi contains cell-independent grid ROI focus scores transferred to cells)
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    """
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))

    # Plot 1: Ranked scatter plot with density coloring
    ax = axes[0, 0]
    # Calculate ranks (higher rank = higher focus score)
    ccfs_ranks = new_df["CCFS_DAPI"].rank(method="average")
    roi_ranks = new_df["DAPI_RFSnorm_roi"].rank(method="average")

    # Create scatter plot with density coloring (subsample for KDE performance)
    try:
        from scipy.stats import gaussian_kde

        valid = ccfs_ranks.dropna().index.intersection(roi_ranks.dropna().index)
        x_vals, y_vals = ccfs_ranks[valid].values, roi_ranks[valid].values
        max_kde_pts = 20000
        if len(x_vals) > max_kde_pts:
            rng = np.random.default_rng(42)
            sidx = rng.choice(len(x_vals), max_kde_pts, replace=False)
            xy_sub = np.vstack([x_vals[sidx], y_vals[sidx]])
            kde = gaussian_kde(xy_sub)
            z = kde(np.vstack([x_vals, y_vals]))
        else:
            z = gaussian_kde(np.vstack([x_vals, y_vals]))(np.vstack([x_vals, y_vals]))
        idx = z.argsort()
        scatter = ax.scatter(
            x_vals[idx],
            y_vals[idx],
            c=z[idx],
            s=1,
            cmap="viridis",
            alpha=0.6,
            rasterized=True,
        )
        plt.colorbar(scatter, ax=ax, label="Density")
    except (ImportError, Exception):
        # Fallback if scipy not available or density calculation fails
        ax.scatter(ccfs_ranks, roi_ranks, alpha=0.3, s=1, marker="o", rasterized=True)

    ax.set_xlabel("CCFS_DAPI Rank (Nuclei-based)", fontsize=12)
    ax.set_ylabel("DAPI_RFSnorm_roi Rank (Cell-Independent Tile-based)", fontsize=12)
    ax.set_title("CCFS vs Cell-Independent Tile Focus Score (Ranked)", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.7)

    # Add correlation coefficient (Spearman rank-based correlation)
    correlation = new_df["CCFS_DAPI"].corr(
        new_df["DAPI_RFSnorm_roi"], method="spearman"
    )
    ax.text(
        0.05,
        0.95,
        f"Spearman Correlation: {correlation:.4f}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # Plot 2: Density plot overlay
    ax = axes[0, 1]
    ax.scatter(
        new_df["CCFS_DAPI"],
        new_df["DAPI_RFSnorm_roi"],
        alpha=0.1,
        s=0.5,
        marker="o",
        c="blue",
        label="Data points",
        rasterized=True,
    )
    # Add 2D density contour (subsample for KDE performance)
    try:
        from scipy.stats import gaussian_kde

        valid_mask = new_df["CCFS_DAPI"].notna() & new_df["DAPI_RFSnorm_roi"].notna()
        x_vals = new_df.loc[valid_mask, "CCFS_DAPI"].values
        y_vals = new_df.loc[valid_mask, "DAPI_RFSnorm_roi"].values
        max_kde_pts = 20000
        if len(x_vals) > max_kde_pts:
            rng = np.random.default_rng(42)
            sidx = rng.choice(len(x_vals), max_kde_pts, replace=False)
            kde = gaussian_kde(np.vstack([x_vals[sidx], y_vals[sidx]]))
            z = kde(np.vstack([x_vals, y_vals]))
        else:
            z = gaussian_kde(np.vstack([x_vals, y_vals]))(np.vstack([x_vals, y_vals]))
        idx = z.argsort()
        ax.scatter(
            x_vals[idx],
            y_vals[idx],
            c=z[idx],
            s=1,
            cmap="viridis",
            alpha=0.6,
            rasterized=True,
        )
    except (ImportError, Exception):
        pass  # Skip density overlay if unavailable
    ax.set_xlabel("CCFS_DAPI (Nuclei-based)", fontsize=12)
    ax.set_ylabel("DAPI_RFSnorm_roi (Cell-Independent Tile-based)", fontsize=12)
    ax.set_title("CCFS vs Cell-Independent Tile Focus Score (Density)", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.7)

    # Plot 3: Agreement/disagreement classification
    ax = axes[1, 0]
    # Prefer GMM 2D if available, otherwise use threshold-based
    has_gmm_2d = "is_blurred_gmm_2d_roi" in new_df.columns
    has_threshold = "is_blurred_roi" in new_df.columns

    if "is_low_nuclear_texture" in new_df.columns and (has_gmm_2d or has_threshold):
        # Use GMM 2D if available, otherwise threshold-based
        if has_gmm_2d:
            roi_blurred_col = "is_blurred_gmm_2d_roi"
            roi_high_focus = ~new_df["is_blurred_gmm_2d_roi"]
            method_label = "2D GMM"
        else:
            roi_blurred_col = "is_blurred_roi"
            roi_high_focus = new_df["is_high_focus_roi"]
            method_label = "Threshold"

        # Create agreement categories
        agreement = new_df["is_low_nuclear_texture"] == new_df[roi_blurred_col]
        agree_high = new_df["is_high_nuclear_texture"] & roi_high_focus
        agree_low = new_df["is_low_nuclear_texture"] & new_df[roi_blurred_col]
        disagree = ~agreement

        # Plot
        ax.scatter(
            new_df.loc[agree_high, "CCFS_DAPI"],
            new_df.loc[agree_high, "DAPI_RFSnorm_roi"],
            alpha=0.3,
            s=1,
            c="green",
            label="Both high focus",
            marker="o",
            rasterized=True,
        )
        ax.scatter(
            new_df.loc[agree_low, "CCFS_DAPI"],
            new_df.loc[agree_low, "DAPI_RFSnorm_roi"],
            alpha=0.3,
            s=1,
            c="red",
            label="Both blurred",
            marker="o",
            rasterized=True,
        )
        ax.scatter(
            new_df.loc[disagree, "CCFS_DAPI"],
            new_df.loc[disagree, "DAPI_RFSnorm_roi"],
            alpha=0.5,
            s=2,
            c="orange",
            label="Disagree",
            marker="x",
            rasterized=True,
        )
        ax.set_xlabel("CCFS_DAPI (Nuclei-based)", fontsize=12)
        ax.set_ylabel("DAPI_RFSnorm_roi (Cell-Independent Tile-based)", fontsize=12)
        ax.set_title(
            f"Classification Agreement (CCFS vs Tile-based {method_label})", fontsize=14
        )
        ax.legend(fontsize=10)
        ax.grid(True, linestyle="--", alpha=0.7)

        # Add agreement percentage
        agreement_rate = agreement.sum() / len(new_df) * 100
        ax.text(
            0.05,
            0.95,
            f"Agreement: {agreement_rate:.2f}%",
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

    # Plot 4: Distribution comparison
    ax = axes[1, 1]
    ax.hist(
        new_df["CCFS_DAPI"].dropna(),
        bins=50,
        alpha=0.5,
        label="CCFS_DAPI (Nuclei-based)",
        color="blue",
        density=True,
    )
    ax.hist(
        new_df["DAPI_RFSnorm_roi"].dropna(),
        bins=50,
        alpha=0.5,
        label="DAPI_RFSnorm_roi (Cell-Independent Tile)",
        color="red",
        density=True,
    )
    ax.set_xlabel("Focus Score", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Distribution Comparison (CCFS vs Cell-Independent Tile)", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.savefig(
        figures_dir / "ccfs_vs_roi_comparison.pdf", dpi=300, bbox_inches="tight"
    )
    plt.savefig(
        figures_dir / "ccfs_vs_roi_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data as CSV (including ranks)
    df_comparison = new_df[["CCFS_DAPI", "DAPI_RFS_roi", "DAPI_RFSnorm_roi"]].copy()
    # Add ranks for analysis
    df_comparison["CCFS_DAPI_rank"] = new_df["CCFS_DAPI"].rank(method="average")
    df_comparison["DAPI_RFSnorm_roi_rank"] = new_df["DAPI_RFSnorm_roi"].rank(
        method="average"
    )

    # Add classification columns - prefer GMM 2D, include threshold-based if available
    if "is_low_nuclear_texture" in new_df.columns:
        df_comparison["is_low_nuclear_texture"] = new_df["is_low_nuclear_texture"]

        # GMM 2D classification (preferred)
        if "is_blurred_gmm_2d_roi" in new_df.columns:
            df_comparison["is_blurred_gmm_2d_roi"] = new_df["is_blurred_gmm_2d_roi"]
            df_comparison["classification_agreement_gmm_2d"] = (
                new_df["is_low_nuclear_texture"] == new_df["is_blurred_gmm_2d_roi"]
            )
            if "blur_prob_gmm_2d_roi" in new_df.columns:
                df_comparison["blur_prob_gmm_2d_roi"] = new_df["blur_prob_gmm_2d_roi"]

        # Threshold-based classification (for comparison)
        if "is_blurred_roi" in new_df.columns:
            df_comparison["is_blurred_roi"] = new_df["is_blurred_roi"]
            df_comparison["classification_agreement"] = (
                new_df["is_low_nuclear_texture"] == new_df["is_blurred_roi"]
            )

    df_comparison.to_csv(figures_source_dir / "ccfs_vs_roi_comparison.csv", index=False)


def plot_spatial_comparison(new_df, myData, figures_dir, figures_source_dir):
    """
    Create spatial comparison plots showing both nuclei-based and tile-based focus scores.

    Parameters:
    -----------
    new_df : pandas DataFrame
        DataFrame containing spatial and focus score data
    myData : pandas DataFrame
        DataFrame with nuclei-based measurements (for spatial coordinates)
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    """
    # Compute figure size from data extents to avoid empty white space
    _x_range = (
        myData["centroid-1"].max() - myData["centroid-1"].min()
        if "centroid-1" in myData.columns
        else 1
    )
    _y_range = (
        myData["centroid-0"].max() - myData["centroid-0"].min()
        if "centroid-0" in myData.columns
        else 1
    )
    _data_aspect = _y_range / _x_range if _x_range > 0 else 1.0
    _pw = 7
    _ph = max(4, _pw * _data_aspect)
    fig, axes = plt.subplots(1, 2, figsize=(2 * _pw + 2, _ph))

    # Plot 1: Nuclei-based CCFS spatial
    ax = axes[0]
    if "centroid-1" in myData.columns and "centroid-0" in myData.columns:
        # Filter out NaN values for plotting
        valid_mask = myData["CCFS_DAPI"].notna()
        if valid_mask.sum() > 0:
            # Use CCFS_DAPI directly (not negated) with auto-scaling
            ccfs_values = myData.loc[valid_mask, "CCFS_DAPI"]
            scatter = ax.scatter(
                myData.loc[valid_mask, "centroid-1"],
                -myData.loc[valid_mask, "centroid-0"],
                s=0.1,
                c=ccfs_values,
                cmap="viridis",
                rasterized=True,
            )
            ax.set_title("Nuclei-based Focus Score (CCFS_DAPI)", fontsize=14)
            ax.set_facecolor("black")
            ax.set_aspect("equal")
            plt.colorbar(scatter, ax=ax, label="CCFS_DAPI")
        else:
            ax.text(
                0.5,
                0.5,
                "No valid CCFS_DAPI data",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=14,
                bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
            )
            ax.set_facecolor("black")
            ax.set_aspect("equal")

    # Plot 2: Tile-based RFS spatial (colored by GMM 2D classification if available)
    ax = axes[1]
    if "DAPI_RFSnorm_roi" in new_df.columns:
        # Prefer GMM 2D classification for coloring, otherwise use raw focus scores
        has_gmm_2d = "is_blurred_gmm_2d_roi" in new_df.columns
        valid_mask = new_df["DAPI_RFSnorm_roi"].notna()

        if has_gmm_2d:
            # Color by GMM 2D classification (blurred vs in-focus)
            gmm_valid_mask = valid_mask & new_df["is_blurred_gmm_2d_roi"].notna()
            if gmm_valid_mask.sum() > 0:
                # Color: red for blurred, blue for in-focus
                colors = [
                    "red" if blurred else "blue"
                    for blurred in new_df.loc[gmm_valid_mask, "is_blurred_gmm_2d_roi"]
                ]
                ax.scatter(
                    new_df.loc[gmm_valid_mask, "x"],
                    -new_df.loc[gmm_valid_mask, "y"],
                    s=0.1,
                    c=colors,
                    alpha=0.5,
                    rasterized=True,
                )
                ax.set_title(
                    "Tile-based Focus Score (GMM 2D Classification)", fontsize=14
                )
                # Add legend
                from matplotlib.patches import Patch

                legend_elements = [
                    Patch(facecolor="blue", label="In-Focus (2D GMM)"),
                    Patch(facecolor="red", label="Blurred (2D GMM)"),
                ]
                ax.legend(handles=legend_elements, fontsize=10, loc="upper right")
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No valid GMM 2D classification data",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=14,
                    bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
                )
        else:
            # Fallback: use raw focus scores with color scale
            if valid_mask.sum() > 0:
                roi_values = new_df.loc[valid_mask, "DAPI_RFSnorm_roi"]
                scatter = ax.scatter(
                    new_df.loc[valid_mask, "x"],
                    -new_df.loc[valid_mask, "y"],
                    s=0.1,
                    c=roi_values,
                    cmap="viridis",
                    rasterized=True,
                )
                ax.set_title("Tile-based Focus Score (DAPI_RFSnorm_roi)", fontsize=14)
                plt.colorbar(scatter, ax=ax, label="DAPI_RFSnorm_roi")
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No valid DAPI_RFSnorm_roi data",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=14,
                    bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.5),
                )

        ax.set_facecolor("black")
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(
        figures_dir / "spatial_comparison_nuclei_vs_roi.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        figures_dir / "spatial_comparison_nuclei_vs_roi.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data as CSV
    if "centroid-1" in myData.columns:
        df_spatial_comparison = pd.DataFrame(
            {
                "x_nuclei": myData["centroid-1"],
                "y_nuclei": -myData["centroid-0"],
                "CCFS_DAPI": myData["CCFS_DAPI"],
            }
        )
        if "DAPI_RFSnorm_roi" in new_df.columns:
            # Try to merge on coordinates
            df_spatial_comparison = df_spatial_comparison.merge(
                new_df[["x", "y", "DAPI_RFSnorm_roi"]],
                left_on=["x_nuclei", "y_nuclei"],
                right_on=["x", "y"],
                how="left",
            )
        df_spatial_comparison.to_csv(
            figures_source_dir / "spatial_comparison_nuclei_vs_roi.csv", index=False
        )


def save_cell_qc_metrics(new_df, outdir, roi_size=None):
    """Save QC metrics to CSV and JSON"""

    # Move cell_id to first position
    cell_id_col = new_df.pop("cell_id")
    new_df.insert(0, "cell_id", cell_id_col)
    # Save full dataset
    new_df.to_csv(outdir / "image_qc_cell_metrics.csv", index=False)

    # Generate QC summary statistics
    qc_metrics = {
        "total_cells": len(new_df),
        "cells_with_transcripts": len(new_df[new_df["transcript_counts"] > 0]),
        "mean_transcript_count": float(new_df["transcript_counts"].mean()),
        "median_transcript_count": float(new_df["transcript_counts"].median()),
        "std_transcript_count": float(new_df["transcript_counts"].std()),
        "mean_ccfs_dapi": float(new_df["CCFS_DAPI"].mean()),
        "cells_high_nuclear_texture": len(new_df[new_df["is_high_nuclear_texture"]]),
        "cells_low_nuclear_texture": len(new_df[new_df["is_low_nuclear_texture"]]),
        "cells_near_edge": len(new_df[new_df["is_near_edge"]]),
        "cells_near_holes": len(new_df[new_df["is_near_hole"]]),
        "cells_with_dense_intensity_regions": len(
            new_df[new_df["has_dense_intensity_regions"]]
        ),
        "unique_dense_intensity_regions": sorted(
            new_df["Dense-Intensity-Region-ID"].unique().tolist()
        )
        if "Dense-Intensity-Region-ID" in new_df.columns
        else [],
        "total_dense_intensity_regions": len(
            new_df[new_df["Dense-Intensity-Region-ID"] > 0][
                "Dense-Intensity-Region-ID"
            ].unique()
        )
        if "Dense-Intensity-Region-ID" in new_df.columns
        else 0,
        "clusters_present": sorted(new_df["Cluster_kmeans10"].unique().tolist()),
        "segmentation_methods": sorted(new_df["segmentation_method"].unique().tolist()),
    }

    # Add tile-based metrics if available
    if "DAPI_RFS_roi" in new_df.columns:
        # Save ROI size used for this analysis
        if roi_size is not None:
            qc_metrics["roi_size_used"] = int(roi_size)
        qc_metrics["mean_dapi_rfs_roi"] = float(new_df["DAPI_RFS_roi"].mean())
        qc_metrics["median_dapi_rfs_roi"] = float(new_df["DAPI_RFS_roi"].median())
        qc_metrics["mean_dapi_rfsnorm_roi"] = float(new_df["DAPI_RFSnorm_roi"].mean())
        qc_metrics["median_dapi_rfsnorm_roi"] = float(
            new_df["DAPI_RFSnorm_roi"].median()
        )

        # GMM 2D metrics (preferred method)
        if "is_blurred_gmm_2d_roi" in new_df.columns:
            qc_metrics["cells_high_focus_gmm_2d_roi"] = len(
                new_df[~new_df["is_blurred_gmm_2d_roi"]]
            )
            qc_metrics["cells_blurred_gmm_2d_roi"] = len(
                new_df[new_df["is_blurred_gmm_2d_roi"]]
            )
            if "blur_prob_gmm_2d_roi" in new_df.columns:
                valid_probs = new_df["blur_prob_gmm_2d_roi"].dropna()
                if len(valid_probs) > 0:
                    qc_metrics["mean_blur_prob_gmm_2d_roi"] = float(valid_probs.mean())
                    qc_metrics["median_blur_prob_gmm_2d_roi"] = float(
                        valid_probs.median()
                    )

        # Threshold-based metrics (for comparison/backward compatibility)
        if "is_blurred_roi" in new_df.columns:
            qc_metrics["cells_high_focus_roi"] = len(
                new_df[new_df["is_high_focus_roi"]]
            )
            qc_metrics["cells_blurred_roi"] = len(new_df[new_df["is_blurred_roi"]])

        # Comparison metrics between nuclei-based and tile-based methods
        # Calculate Spearman rank-based correlation between CCFS_DAPI and DAPI_RFSnorm_roi
        # Spearman is better suited for different scales and non-linear relationships
        if "CCFS_DAPI" in new_df.columns and "DAPI_RFSnorm_roi" in new_df.columns:
            correlation = new_df["CCFS_DAPI"].corr(
                new_df["DAPI_RFSnorm_roi"], method="spearman"
            )
            qc_metrics["ccfs_vs_rfs_correlation"] = (
                float(correlation) if not np.isnan(correlation) else None
            )
            qc_metrics["ccfs_vs_rfs_correlation_method"] = "spearman"

            # Calculate agreement rate - prefer GMM 2D, include threshold-based for comparison
            if "is_low_nuclear_texture" in new_df.columns:
                # GMM 2D agreement (preferred)
                if "is_blurred_gmm_2d_roi" in new_df.columns:
                    agreement_gmm_2d = (
                        new_df["is_low_nuclear_texture"]
                        == new_df["is_blurred_gmm_2d_roi"]
                    ).sum()
                    agreement_rate_gmm_2d = agreement_gmm_2d / len(new_df)
                    qc_metrics["classification_agreement_gmm_2d"] = float(
                        agreement_rate_gmm_2d
                    )
                    qc_metrics["classification_agreement_gmm_2d_count"] = int(
                        agreement_gmm_2d
                    )
                    qc_metrics["classification_disagreement_gmm_2d_count"] = int(
                        len(new_df) - agreement_gmm_2d
                    )

                # Threshold-based agreement (for comparison)
                if "is_blurred_roi" in new_df.columns:
                    agreement = (
                        new_df["is_low_nuclear_texture"] == new_df["is_blurred_roi"]
                    ).sum()
                    agreement_rate = agreement / len(new_df)
                    qc_metrics["classification_agreement"] = float(agreement_rate)
                    qc_metrics["classification_agreement_count"] = int(agreement)
                    qc_metrics["classification_disagreement_count"] = int(
                        len(new_df) - agreement
                    )

    # Save metrics
    with open(outdir / "image_qc_cell_metrics.json", "w") as f:
        json.dump(qc_metrics, f, indent=2)

    # Save ROI size to simple text file for easy retrieval
    if roi_size is not None:
        with open(outdir / "roi_size.txt", "w") as f:
            f.write(f"{roi_size}\n")

    logging.info("\n=== IMAGE QC SUMMARY ===")
    logging.info(f"Total cells analyzed: {qc_metrics['total_cells']:,}")
    logging.info(f"Cells with transcripts: {qc_metrics['cells_with_transcripts']:,}")
    logging.info(f"Mean transcript count: {qc_metrics['mean_transcript_count']:.1f}")
    logging.info(f"Mean CCFS DAPI: {qc_metrics['mean_ccfs_dapi']:.6f}")
    logging.info(
        f"Cells high nuclear texture: {qc_metrics['cells_high_nuclear_texture']:,}"
    )
    logging.info(
        f"Cells low nuclear texture: {qc_metrics['cells_low_nuclear_texture']:,}"
    )
    logging.info(f"Cells near edge: {qc_metrics['cells_near_edge']:,}")
    logging.info(f"Cells near holes: {qc_metrics['cells_near_holes']:,}")
    logging.info(
        f"Cells with dense intensity regions: {qc_metrics['cells_with_dense_intensity_regions']:,}"
    )
    if (
        "total_dense_intensity_regions" in qc_metrics
        and qc_metrics["total_dense_intensity_regions"] > 0
    ):
        logging.info(
            f"Total unique dense intensity regions detected: {qc_metrics['total_dense_intensity_regions']}"
        )
        logging.info(
            f"  (Region IDs: {qc_metrics['unique_dense_intensity_regions'][:10]}{'...' if len(qc_metrics['unique_dense_intensity_regions']) > 10 else ''})"
        )
    logging.info(f"Clusters present: {qc_metrics['clusters_present']}")
    logging.info(f"Segmentation methods: {qc_metrics['segmentation_methods']}")

    # Print tile-based metrics if available
    if "mean_dapi_rfs_roi" in qc_metrics:
        logging.info("\n=== TILE-BASED FOCUS SCORE SUMMARY ===")
        logging.info(f"Mean DAPI RFS (tile): {qc_metrics['mean_dapi_rfs_roi']:.6f}")
        logging.info(
            f"Mean DAPI RFS normalized (tile): {qc_metrics['mean_dapi_rfsnorm_roi']:.6f}"
        )

        # GMM 2D metrics (preferred)
        if "cells_blurred_gmm_2d_roi" in qc_metrics:
            logging.info("\n--- 2D GMM Classification (Preferred) ---")
            logging.info(
                f"Cells with high focus (2D GMM): {qc_metrics['cells_high_focus_gmm_2d_roi']:,}"
            )
            logging.info(
                f"Cells blurred (2D GMM): {qc_metrics['cells_blurred_gmm_2d_roi']:,}"
            )
            if "mean_blur_prob_gmm_2d_roi" in qc_metrics:
                logging.info(
                    f"Mean blur probability (2D GMM): {qc_metrics['mean_blur_prob_gmm_2d_roi']:.4f}"
                )

        # Threshold-based metrics (for comparison)
        if "cells_blurred_roi" in qc_metrics:
            logging.info("\n--- Threshold-based Classification (Comparison) ---")
            logging.info(
                f"Cells with high focus (Threshold): {qc_metrics['cells_high_focus_roi']:,}"
            )
            logging.info(
                f"Cells blurred (Threshold): {qc_metrics['cells_blurred_roi']:,}"
            )

        if (
            "ccfs_vs_rfs_correlation" in qc_metrics
            and qc_metrics["ccfs_vs_rfs_correlation"] is not None
        ):
            logging.info("\n=== METHOD COMPARISON ===")
            logging.info(
                f"Spearman Correlation (CCFS vs RFS): {qc_metrics['ccfs_vs_rfs_correlation']:.4f}"
            )

            # GMM 2D agreement (preferred)
            if "classification_agreement_gmm_2d" in qc_metrics:
                logging.info("\n--- CCFS vs 2D GMM Agreement (Preferred) ---")
                logging.info(
                    f"Classification agreement: {qc_metrics['classification_agreement_gmm_2d'] * 100:.2f}%"
                )
                logging.info(
                    f"  - Agreeing cells: {qc_metrics['classification_agreement_gmm_2d_count']:,}"
                )
                logging.info(
                    f"  - Disagreeing cells: {qc_metrics['classification_disagreement_gmm_2d_count']:,}"
                )

            # Threshold-based agreement (for comparison)
            if "classification_agreement" in qc_metrics:
                logging.info("\n--- CCFS vs Threshold Agreement (Comparison) ---")
                logging.info(
                    f"Classification agreement: {qc_metrics['classification_agreement'] * 100:.2f}%"
                )
                logging.info(
                    f"  - Agreeing cells: {qc_metrics['classification_agreement_count']:,}"
                )
                logging.info(
                    f"  - Disagreeing cells: {qc_metrics['classification_disagreement_count']:,}"
                )


def generate_cell_figures(
    data,
    new_df,
    myData,
    figures_dir,
    figures_source_dir,
    roi_threshold=None,
    roi_intensity_threshold=None,
    ccfs_low_texture_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
):
    """
    Generate cell-based figures from merged data using multithreading.

    Parameters:
    -----------
    data : dict
        Dictionary with paths and directories
    new_df : pandas DataFrame
        Merged DataFrame with cell data and ROI mappings
    myData : pandas DataFrame
        Nuclei-based measurements DataFrame
    figures_dir : Path
        Directory to save figures
    figures_source_dir : Path
        Directory to save source data
    roi_threshold : float, optional
        Tile-based blur threshold for display
    roi_intensity_threshold : float, optional
        Tile intensity threshold for display
    """
    figures_source_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Generating cell-based figures...")

    def _plot1_nuclear_texture_proportions():
        logging.info("  - Nuclear texture proportions (CCFS)...")
        plot_nuclear_texture_proportions(
            new_df,
            figures_dir,
            figures_source_dir,
            texture_threshold=ccfs_low_texture_threshold,
        )

    def _plot2_blur_proportions_roi():
        if (
            "is_blurred_gmm_2d_roi" in new_df.columns
            or "is_blurred_roi" in new_df.columns
        ):
            logging.info("  - Blur score proportions (tile-based)...")
            roi_threshold_display = roi_threshold if roi_threshold is not None else -1.0
            plot_tile_blur_proportions_roi._intensity_threshold = (
                roi_intensity_threshold if roi_intensity_threshold is not None else 20.0
            )
            plot_tile_blur_proportions_roi(
                new_df,
                figures_dir,
                figures_source_dir,
                roi_threshold=roi_threshold_display,
            )

    def _plot3_nuclear_texture_density():
        logging.info("  - Nuclear texture density...")
        plot_nuclear_texture_density(
            new_df,
            figures_dir,
            figures_source_dir,
            ccfs_low_texture_threshold=ccfs_low_texture_threshold,
        )

    def _plot5_ccfs_vs_roi():
        if "DAPI_RFSnorm_roi" in new_df.columns:
            logging.info("  - CCFS vs tile comparison...")
            plot_ccfs_vs_roi_comparison(new_df, figures_dir, figures_source_dir)

    def _plot6_spatial_comparison():
        logging.info("  - Spatial comparison...")
        plot_spatial_comparison(new_df, myData, figures_dir, figures_source_dir)

    def _plot7_cell_focus_distribution():
        logging.info("  - Cell-level focus score distribution...")
        plot_cell_focus_distribution(new_df, figures_dir, figures_source_dir)

    def _plot8_gmm_focus_vs_transcripts():
        logging.info("  - GMM focus vs transcripts...")
        plot_gmm_focus_vs_transcripts(new_df, figures_dir, figures_source_dir)

    tasks = [
        _plot1_nuclear_texture_proportions,
        _plot2_blur_proportions_roi,
        _plot3_nuclear_texture_density,
        _plot5_ccfs_vs_roi,
        _plot6_spatial_comparison,
        _plot7_cell_focus_distribution,
        _plot8_gmm_focus_vs_transcripts,
    ]

    # Use multiprocessing.Process with fork context -- matplotlib is NOT thread-safe.
    # fork is intentional: this code runs exclusively in Linux/Docker containers where
    # fork is safe, and the closures capture large numpy arrays that cannot be pickled
    # (required by spawn/forkserver). Do NOT change to spawn/forkserver without
    # restructuring the closure-based task definitions.
    max_workers = min(len(tasks), os.cpu_count() or 4)
    mp_ctx = multiprocessing.get_context("fork")

    # Run in batches to limit concurrent processes
    failed = []
    for batch_start in range(0, len(tasks), max_workers):
        batch = tasks[batch_start : batch_start + max_workers]
        processes = []
        for fn in batch:
            p = mp_ctx.Process(target=_run_figure_task, args=(fn,))
            p.start()
            processes.append((p, fn.__name__))
        for p, name in processes:
            p.join()
            if p.exitcode != 0:
                failed.append(name)

    if failed:
        raise RuntimeError(f"Cell figure generation failed for: {', '.join(failed)}")

    logging.info("All cell-based figures generated successfully!")


# ===== QUARTO FIGURE GENERATION (from image_qc_processing.py) =====


def generate_all_figures(
    data,
    df_spatial,
    new_df,
    myData,
    small0,
    small1,
    small2,
    distance_map,
    distance_map2,
    whole_sample,
    holes,
    artefacts,
    ccfs_low_texture_threshold=DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD,
):
    """Generate ALL 13 Quarto-required figures using multithreading, and save the data used for each plot as CSV."""
    figures_dir = data["figures_dir"]
    figures_source_dir = figures_dir / "figures_source"
    figures_source_dir.mkdir(parents=True, exist_ok=True)
    # if df_spatial['In-Area-with-Artefact'] have a single unique value, use that value, otherwise use 0.5
    if len(df_spatial["In-Area-with-Artefact"].unique()) == 1:
        iawa_vmax = df_spatial["In-Area-with-Artefact"].unique()[0]
    else:
        iawa_vmax = 0.5

    def _fig1_distance_edge():
        logging.info("Generating Figure 1: Distance map (edge)...")
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=True, sharey=True)
        ax.imshow(distance_map, cmap="gray")
        plt.tight_layout()
        plt.savefig(figures_dir / "distance_map_edge.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(distance_map).to_csv(
            figures_source_dir / "distance_map_edge.csv", index=False, header=False
        )

    def _fig2_distance_holes():
        logging.info("Generating Figure 2: Distance map (holes)...")
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=True, sharey=True)
        ax.imshow(distance_map2, cmap="gray")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "distance_map_holes.pdf", dpi=300, bbox_inches="tight"
        )
        plt.savefig(
            figures_dir / "distance_map_holes.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        pd.DataFrame(distance_map2).to_csv(
            figures_source_dir / "distance_map_holes.csv", index=False, header=False
        )

    def _fig3_morphology_overview():
        logging.info("Generating Figure 3: Morphology overview...")
        r = (small0.shape[0] / small0.shape[1]) * 12
        fig, ax = plt.subplots(2, 2, figsize=(12, r))
        ax[0, 0].imshow(small0, cmap="Greys_r", vmax=np.percentile(small0, 99))
        ax[0, 0].set_title("DAPI")
        ax[0, 0].set_aspect("equal")
        if small1 is not None:
            ax[0, 1].imshow(small1, cmap="Greys_r", vmax=np.percentile(small1, 99))
            ax[0, 1].set_title("Boundary")
            ax[0, 1].set_aspect("equal")
        else:
            ax[0, 1].set_title("Boundary (not available)")
            ax[0, 1].axis("off")
        if small2 is not None:
            ax[1, 0].imshow(small2, cmap="Greys_r", vmax=np.percentile(small2, 99))
            ax[1, 0].set_title("Interior")
            ax[1, 0].set_aspect("equal")
        else:
            ax[1, 0].set_title("Interior (not available)")
            ax[1, 0].axis("off")
        ax[1, 1].set_title("Sample artifacts map")
        ax[1, 1].imshow(color.label2rgb(artefacts, bg_label=0))
        ax[1, 1].set_aspect("equal")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "morphology_overview.pdf", dpi=300, bbox_inches="tight"
        )
        plt.savefig(
            figures_dir / "morphology_overview.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        pd.DataFrame(small0).to_csv(
            figures_source_dir / "morphology_overview_DAPI.csv",
            index=False,
            header=False,
        )
        if small1 is not None:
            pd.DataFrame(small1).to_csv(
                figures_source_dir / "morphology_overview_Boundary.csv",
                index=False,
                header=False,
            )
        if small2 is not None:
            pd.DataFrame(small2).to_csv(
                figures_source_dir / "morphology_overview_Interior.csv",
                index=False,
                header=False,
            )
        pd.DataFrame(artefacts).to_csv(
            figures_source_dir / "morphology_overview_Artefacts.csv",
            index=False,
            header=False,
        )

    def _fig4_sample_qc_metrics():
        logging.info("Generating Figure 4: Sample QC metrics...")
        idx4 = _subsample_idx(np.ones(len(df_spatial), dtype=bool))
        fig, ax = plt.subplots(2, 2, figsize=(12, 10))
        ax[0, 0].imshow(small0, cmap="Greys_r", vmax=np.percentile(small0, 99))
        ax[0, 0].set_title("DAPI")
        ax[0, 0].set_aspect("equal")
        ax[0, 1].scatter(
            df_spatial["x"].iloc[idx4],
            -df_spatial["y"].iloc[idx4],
            c=df_spatial["Distance-to-edge"].iloc[idx4],
            s=0.1,
            marker="o",
            rasterized=True,
        )
        ax[0, 1].set_title("Distance to sample Edge")
        ax[0, 1].set_aspect("equal")
        ax[0, 1].set_facecolor("black")
        ax[1, 0].scatter(
            df_spatial["x"].iloc[idx4],
            -df_spatial["y"].iloc[idx4],
            c=df_spatial["Distance-to-nearest-hole"].iloc[idx4],
            s=0.1,
            marker="o",
            rasterized=True,
        )
        ax[1, 0].set_title("Distance to nearest hole")
        ax[1, 0].set_aspect("equal")
        ax[1, 0].set_facecolor("black")
        ax[1, 1].scatter(
            df_spatial["x"].iloc[idx4],
            -df_spatial["y"].iloc[idx4],
            c=df_spatial["In-Area-with-Artefact"].iloc[idx4],
            cmap="Blues_r",
            s=0.1,
            marker="o",
            vmax=iawa_vmax,
            rasterized=True,
        )
        ax[1, 1].set_title("Overlap with sample artefact")
        ax[1, 1].set_aspect("equal")
        ax[1, 1].set_facecolor("black")
        plt.tight_layout()
        plt.savefig(figures_dir / "sample_qc_metrics.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        df_sample_qc_metrics = pd.DataFrame(
            {
                "x": df_spatial["x"],
                "y": df_spatial["y"],
                "Distance-to-edge": df_spatial["Distance-to-edge"],
                "Distance-to-nearest-hole": df_spatial["Distance-to-nearest-hole"],
                "In-Area-with-Artefact": df_spatial["In-Area-with-Artefact"],
            }
        )
        df_sample_qc_metrics.to_csv(
            figures_source_dir / "sample_qc_metrics.csv", index=False
        )

    def _fig5_imageqc_masks():
        logging.info("Generating Figure 5: ImageQC masks...")
        fig, ax = plt.subplots(2, 2, figsize=(12, 10))
        ax[0, 0].set_title("DAPI morphology image")
        ax[0, 0].imshow(small0, cmap="Greys_r", vmax=np.percentile(small0, 99))
        ax[0, 0].axis("off")
        ax[0, 1].set_title("Whole Sample mask")
        ax[0, 1].imshow(color.label2rgb(whole_sample, bg_label=0))
        ax[0, 1].axis("off")
        ax[1, 0].set_title("Holes in sample")
        ax[1, 0].imshow(color.label2rgb(holes, bg_label=0))
        ax[1, 0].axis("off")
        ax[1, 1].set_title("Sample artefacts")
        ax[1, 1].imshow(color.label2rgb(artefacts, bg_label=0))
        ax[1, 1].axis("off")
        plt.tight_layout()
        plt.savefig(figures_dir / "imageqc_masks.pdf", dpi=300, bbox_inches="tight")
        plt.savefig(figures_dir / "imageqc_masks.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(small0).to_csv(
            figures_source_dir / "imageqc_masks_DAPI.csv", index=False, header=False
        )
        pd.DataFrame(whole_sample).to_csv(
            figures_source_dir / "imageqc_masks_WholeSample.csv",
            index=False,
            header=False,
        )
        pd.DataFrame(holes).to_csv(
            figures_source_dir / "imageqc_masks_Holes.csv", index=False, header=False
        )
        pd.DataFrame(artefacts).to_csv(
            figures_source_dir / "imageqc_masks_Artefacts.csv",
            index=False,
            header=False,
        )

    def _fig6_ccfs_spatial():
        logging.info("Generating Figure 6: CCFS Spatial...")
        idx6 = _subsample_idx(np.ones(len(myData), dtype=bool))
        c_vals = -myData["CCFS_DAPI"].iloc[idx6]
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.scatter(
            myData["centroid-1"].iloc[idx6],
            -myData["centroid-0"].iloc[idx6],
            s=0.1,
            c=c_vals,
            cmap="viridis",
            vmin=c_vals.min(),
            vmax=max(c_vals.max(), c_vals.min() + 1e-6),
            rasterized=True,
        )
        ax.set_title("Calculated Cell Focus Score")
        ax.set_facecolor("black")
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.savefig(figures_dir / "ccfs_spatial.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        df_ccfs_spatial = pd.DataFrame(
            {
                "centroid-1": myData["centroid-1"],
                "centroid-0": myData["centroid-0"],
                "CCFS_DAPI": myData["CCFS_DAPI"],
            }
        )
        df_ccfs_spatial.to_csv(figures_source_dir / "ccfs_spatial.csv", index=False)

    def _fig7_ccfs_thresholded():
        logging.info("Generating Figure 7: CCFS Thresholded...")
        highC = myData[myData["is_high_nuclear_texture"]]
        lowC = myData[myData["is_low_nuclear_texture"]]
        idx7h = _subsample_idx(np.ones(len(highC), dtype=bool))
        idx7l = _subsample_idx(np.ones(len(lowC), dtype=bool))
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.scatter(
            highC["centroid-1"].iloc[idx7h],
            -highC["centroid-0"].iloc[idx7h],
            s=0.1,
            color="#1B2631",
            rasterized=True,
        )
        ax.scatter(
            lowC["centroid-1"].iloc[idx7l],
            -lowC["centroid-0"].iloc[idx7l],
            s=0.1,
            color="red",
            rasterized=True,
        )
        ax.set_title("Thresholded CCFS (low nuclear texture cells are red!)")
        ax.set_facecolor("black")
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.savefig(figures_dir / "ccfs_thresholded.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        df_ccfs_thresholded = pd.DataFrame(
            {
                "centroid-1": myData["centroid-1"],
                "centroid-0": myData["centroid-0"],
                "CCFS_DAPI": myData["CCFS_DAPI"],
                "is_low_nuclear_texture": myData["is_low_nuclear_texture"],
            }
        )
        df_ccfs_thresholded.to_csv(
            figures_source_dir / "ccfs_thresholded.csv", index=False
        )

    def _fig8_umap_multiple_metrics():
        logging.info("Generating Figure 8: UMAP by multiple metrics...")

        def _safe_scatter(axis, x, y, color_values, title, **kwargs):
            """Scatter plot with guard against constant color data.

            matplotlib Normalize raises ValueError when vmin >= vmax (all
            identical values).  Pad vmax by 1 in that case.  String color
            values (hex codes) are passed through directly.
            """
            c = np.asarray(color_values)
            if c.dtype.kind in ("U", "S", "O"):
                # String colors (e.g. hex "#FFABC3") — no normalization needed
                axis.scatter(
                    x, y, s=0.1, c=color_values, marker="o", rasterized=True, **kwargs
                )
            else:
                c = c.astype(float)
                cmin, cmax = np.nanmin(c), np.nanmax(c)
                if "vmin" not in kwargs:
                    kwargs["vmin"] = cmin
                if "vmax" not in kwargs:
                    kwargs["vmax"] = cmax
                if kwargs["vmin"] >= kwargs["vmax"]:
                    kwargs["vmax"] = kwargs["vmin"] + 1
                axis.scatter(x, y, s=0.1, c=c, marker="o", rasterized=True, **kwargs)
            axis.set_title(title)
            axis.set_facecolor("black")
            axis.set_aspect("equal")

        fig, ax = plt.subplots(2, 2, figsize=(15, 15))
        _safe_scatter(
            ax[0, 0],
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            -new_df["CCFS_DAPI"],
            "UMAP by Nuclear Texture Score",
            cmap="viridis",
            vmin=-0.012,
        )
        _safe_scatter(
            ax[0, 1],
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            new_df["Cluster_kmeans10"],
            "UMAP by Cluster Allocation",
        )
        _safe_scatter(
            ax[1, 0],
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            new_df["segPal"],
            "UMAP by Segmentation method",
        )
        _safe_scatter(
            ax[1, 1],
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            new_df["transcript_counts"],
            "UMAP by Transcript counts",
        )
        plt.tight_layout()
        plt.savefig(
            figures_dir / "umap_multiple_metrics.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        df_umap_multiple_metrics = new_df[
            [
                "UMAP-1",
                "UMAP-2",
                "CCFS_DAPI",
                "Cluster_kmeans10",
                "segPal",
                "transcript_counts",
            ]
        ]
        df_umap_multiple_metrics.to_csv(
            figures_source_dir / "umap_multiple_metrics.csv", index=False
        )

    def _fig9_umap_distance_metrics():
        logging.info("Generating Figure 9: UMAP by distance metrics...")
        highC_new = new_df[new_df["is_high_nuclear_texture"]]
        lowC_new = new_df[new_df["is_low_nuclear_texture"]]
        fig, ax = plt.subplots(2, 2, figsize=(15, 15))
        ax[0, 0].scatter(
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            s=0.1,
            c=-new_df["Distance-to-edge"],
            marker="o",
            rasterized=True,
        )
        ax[0, 0].set_title("UMAP by Distance to edge")
        ax[0, 0].set_facecolor("black")
        ax[0, 0].set_aspect("equal")
        ax[0, 1].scatter(
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            s=0.1,
            c=-new_df["Distance-to-nearest-hole"],
            marker="o",
            rasterized=True,
        )
        ax[0, 1].set_title("UMAP by Distance to nearest hole")
        ax[0, 1].set_facecolor("black")
        ax[0, 1].set_aspect("equal")
        ax[1, 0].scatter(
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            s=0.1,
            c=new_df["In-Area-with-Artefact"],
            cmap="Blues_r",
            marker="o",
            vmax=iawa_vmax,
            rasterized=True,
        )
        ax[1, 0].set_title("UMAP by Overlap with artefact")
        ax[1, 0].set_facecolor("black")
        ax[1, 0].set_aspect("equal")
        ax[1, 1].scatter(
            highC_new["UMAP-1"],
            highC_new["UMAP-2"],
            s=0.1,
            color="#1B2631",
            rasterized=True,
        )
        ax[1, 1].scatter(
            lowC_new["UMAP-1"], lowC_new["UMAP-2"], s=0.1, color="red", rasterized=True
        )
        ax[1, 1].set_title("UMAP by Blurred cells (in red)")
        ax[1, 1].set_facecolor("black")
        ax[1, 1].set_aspect("equal")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "umap_distance_metrics.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        df_umap_distance_metrics = new_df[
            [
                "UMAP-1",
                "UMAP-2",
                "Distance-to-edge",
                "Distance-to-nearest-hole",
                "In-Area-with-Artefact",
                "CCFS_DAPI",
                "is_low_nuclear_texture",
            ]
        ]
        df_umap_distance_metrics.to_csv(
            figures_source_dir / "umap_distance_metrics.csv", index=False
        )

    def _fig10_umap_thresholded_metrics():
        logging.info("Generating Figure 10: UMAP thresholded metrics...")
        highE = new_df[new_df["is_far_from_edge"]]
        lowE = new_df[new_df["is_near_edge"]]
        highH = new_df[new_df["is_far_from_hole"]]
        lowH = new_df[new_df["is_near_hole"]]
        fig, ax = plt.subplots(2, 2, figsize=(15, 15))
        ax[0, 0].scatter(
            highE["UMAP-1"], highE["UMAP-2"], s=0.1, color="#1B2631", rasterized=True
        )
        ax[0, 0].scatter(
            lowE["UMAP-1"], lowE["UMAP-2"], s=0.1, color="red", rasterized=True
        )
        ax[0, 0].set_title("UMAP by Distance to edge")
        ax[0, 0].set_facecolor("black")
        ax[0, 0].set_aspect("equal")
        ax[0, 1].scatter(
            highH["UMAP-1"], highH["UMAP-2"], s=0.1, color="#1B2631", rasterized=True
        )
        ax[0, 1].scatter(
            lowH["UMAP-1"], lowH["UMAP-2"], s=0.1, color="red", rasterized=True
        )
        ax[0, 1].set_title("UMAP by Distance to nearest hole")
        ax[0, 1].set_facecolor("black")
        ax[0, 1].set_aspect("equal")
        ax[1, 0].scatter(
            new_df["UMAP-1"],
            new_df["UMAP-2"],
            s=0.1,
            c=new_df["In-Area-with-Artefact"],
            cmap="Blues_r",
            marker="o",
            vmax=iawa_vmax,
            rasterized=True,
        )
        ax[1, 0].set_title("UMAP by Overlap with artefact")
        ax[1, 0].set_facecolor("black")
        ax[1, 0].set_aspect("equal")
        ax[1, 1].scatter(
            new_df[new_df["is_high_nuclear_texture"]]["UMAP-1"],
            new_df[new_df["is_high_nuclear_texture"]]["UMAP-2"],
            s=0.1,
            color="#1B2631",
            rasterized=True,
        )
        ax[1, 1].scatter(
            new_df[new_df["is_low_nuclear_texture"]]["UMAP-1"],
            new_df[new_df["is_low_nuclear_texture"]]["UMAP-2"],
            s=0.1,
            color="red",
            rasterized=True,
        )
        ax[1, 1].set_title("UMAP by Blurred cells (in red)")
        ax[1, 1].set_facecolor("black")
        ax[1, 1].set_aspect("equal")
        plt.tight_layout()
        plt.savefig(
            figures_dir / "umap_thresholded_metrics.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)
        # Save data as CSV
        df_umap_thresholded_metrics = new_df[
            [
                "UMAP-1",
                "UMAP-2",
                "Distance-to-edge",
                "Distance-to-nearest-hole",
                "In-Area-with-Artefact",
                "CCFS_DAPI",
                "is_near_edge",
                "is_near_hole",
                "is_low_nuclear_texture",
            ]
        ]
        df_umap_thresholded_metrics.to_csv(
            figures_source_dir / "umap_thresholded_metrics.csv", index=False
        )

    def _fig11_nuclear_texture_proportions():
        logging.info("Generating Figure 11: Nuclear Texture Proportions by Cluster...")
        plot_nuclear_texture_proportions(
            new_df,
            figures_dir,
            figures_source_dir,
            texture_threshold=ccfs_low_texture_threshold,
            GROUP_BY_COLUMN="Cluster_kmeans10",
        )

    def _fig12_nuclear_texture_density():
        logging.info("Generating Figure 12: Nuclear Texture Density by Cluster...")
        plot_nuclear_texture_density(
            new_df,
            figures_dir,
            figures_source_dir,
            GROUP_BY_COLUMN="Cluster_kmeans10",
            ccfs_low_texture_threshold=ccfs_low_texture_threshold,
        )

    def _fig13_nuclear_texture_vs_transcripts():
        logging.info(
            "Generating Figure 13: Nuclear Texture vs Transcripts (Log Scale)..."
        )
        plot_nuclear_texture_vs_transcripts(
            new_df,
            figures_dir,
            figures_source_dir,
            log_scale=True,
            ccfs_low_texture_threshold=ccfs_low_texture_threshold,
        )

    def _fig14_cell_focus_distribution():
        logging.info("Generating Figure 14: Cell-level focus score distribution...")
        plot_cell_focus_distribution(new_df, figures_dir, figures_source_dir)

    def _fig15_gmm_blur_proportions_by_cluster():
        if (
            "is_blurred_gmm_2d_roi" in new_df.columns
            or "is_blurred_roi" in new_df.columns
        ):
            logging.info("Generating Figure 15: GMM Blur Proportions by Cluster...")
            plot_tile_blur_proportions_roi(
                new_df,
                figures_dir,
                figures_source_dir,
                GROUP_BY_COLUMN="Cluster_kmeans10",
            )

    def _fig16_gmm_focus_vs_transcripts():
        logging.info("Generating Figure 16: GMM Focus vs Transcripts...")
        plot_gmm_focus_vs_transcripts(new_df, figures_dir, figures_source_dir)

    tasks = [
        _fig1_distance_edge,
        _fig2_distance_holes,
        _fig3_morphology_overview,
        _fig4_sample_qc_metrics,
        _fig5_imageqc_masks,
        _fig6_ccfs_spatial,
        _fig7_ccfs_thresholded,
        _fig8_umap_multiple_metrics,
        _fig9_umap_distance_metrics,
        _fig10_umap_thresholded_metrics,
        _fig11_nuclear_texture_proportions,
        _fig12_nuclear_texture_density,
        _fig13_nuclear_texture_vs_transcripts,
        _fig14_cell_focus_distribution,
        _fig15_gmm_blur_proportions_by_cluster,
        _fig16_gmm_focus_vs_transcripts,
    ]

    # Use multiprocessing.Process with fork context -- matplotlib is NOT thread-safe.
    # fork is intentional: this code runs exclusively in Linux/Docker containers where
    # fork is safe, and the closures capture large numpy arrays that cannot be pickled
    # (required by spawn/forkserver). Do NOT change to spawn/forkserver without
    # restructuring the closure-based task definitions.
    max_workers = min(len(tasks), os.cpu_count() or 4)
    mp_ctx = multiprocessing.get_context("fork")

    # Run in batches to limit concurrent processes
    failed = []
    for batch_start in range(0, len(tasks), max_workers):
        batch = tasks[batch_start : batch_start + max_workers]
        processes = []
        for fn in batch:
            p = mp_ctx.Process(target=_run_figure_task, args=(fn,))
            p.start()
            processes.append((p, fn.__name__))
        for p, name in processes:
            p.join()
            if p.exitcode != 0:
                failed.append(name)

    if failed:
        raise RuntimeError(f"Figure generation failed for: {', '.join(failed)}")

    logging.info("All figures generated successfully!")
    logging.info(f"[OK] Saved figures to {figures_dir}")
    logging.info(f"[OK] Saved source data to {figures_source_dir}")


# ===== SIMPLE QC METRICS (from image_qc_processing.py) =====


def save_simple_qc_metrics(new_df, outdir):
    """Save QC metrics to CSV and JSON"""

    # Move cell_id to first position
    cell_id_col = new_df.pop("cell_id")
    new_df.insert(0, "cell_id", cell_id_col)
    # Save full dataset
    new_df.to_csv(outdir / "image_qc_metrics.csv", index=False)

    # Generate QC summary statistics
    qc_metrics = {
        "total_cells": len(new_df),
        "cells_with_transcripts": len(new_df[new_df["transcript_counts"] > 0]),
        "mean_transcript_count": float(new_df["transcript_counts"].mean()),
        "median_transcript_count": float(new_df["transcript_counts"].median()),
        "std_transcript_count": float(new_df["transcript_counts"].std()),
        "mean_ccfs_dapi": float(new_df["CCFS_DAPI"].mean()),
        "cells_high_nuclear_texture": len(new_df[new_df["is_high_nuclear_texture"]]),
        "cells_low_nuclear_texture": len(new_df[new_df["is_low_nuclear_texture"]]),
        "cells_near_edge": len(new_df[new_df["is_near_edge"]]),
        "cells_near_holes": len(new_df[new_df["is_near_hole"]]),
        "cells_with_artifacts": len(new_df[new_df["has_artifacts"]]),
        "clusters_present": sorted(new_df["Cluster_kmeans10"].unique().tolist()),
        "segmentation_methods": sorted(new_df["segmentation_method"].unique().tolist()),
    }

    # Save metrics
    with open(outdir / "image_qc_metrics.json", "w") as f:
        json.dump(qc_metrics, f, indent=2)

    logging.info("\n=== IMAGE QC SUMMARY ===")
    logging.info(f"Total cells analyzed: {qc_metrics['total_cells']:,}")
    logging.info(f"Cells with transcripts: {qc_metrics['cells_with_transcripts']:,}")
    logging.info(f"Mean transcript count: {qc_metrics['mean_transcript_count']:.1f}")
    logging.info(f"Mean CCFS DAPI: {qc_metrics['mean_ccfs_dapi']:.6f}")
    logging.info(
        f"Cells high nuclear texture: {qc_metrics['cells_high_nuclear_texture']:,}"
    )
    logging.info(
        f"Cells low nuclear texture: {qc_metrics['cells_low_nuclear_texture']:,}"
    )
    logging.info(f"Cells near edge: {qc_metrics['cells_near_edge']:,}")
    logging.info(f"Cells near holes: {qc_metrics['cells_near_holes']:,}")
    logging.info(f"Cells with artifacts: {qc_metrics['cells_with_artifacts']:,}")
    logging.info(f"Clusters present: {qc_metrics['clusters_present']}")
    logging.info(f"Segmentation methods: {qc_metrics['segmentation_methods']}")


# ===== NEW: PIXEL-MAP AGGREGATION FOR CCFS =====


def calculate_ccfs_from_focus_maps(
    focus_maps, cell_masks_zarr, cellseg_mask, xoa_morphology_files
):
    """
    Derive per-cell CCFS from pixel focus maps using scipy.ndimage.mean.

    This replaces the old regionprops-based calculate_ccfs_measurements() for
    DAPI focus scores, but still uses regionprops for boundary/RNA per-cell
    intensities (different channels).

    Args:
        focus_maps: dict from compute_all_focus_maps() with at least 'dapi_focus_map' and 'dapi_mean_map'
        cell_masks_zarr: zarr group with 'masks/0' (nuclear mask) and 'masks/1' (cell mask)
        cellseg_mask: numpy array of cell segmentation mask (from masks/1)
        xoa_morphology_files: list of morphology file paths (for boundary/RNA channels)

    Returns:
        pandas DataFrame with columns: CellID, centroid-0, centroid-1,
        CCFS_DAPI, mean_intensity, area_nucleus, area_cell,
        mean_intensity_Boundary, mean_intensity_IntRNA
    """
    # Load nuclear segmentation mask
    nuclear_mask = np.array(cell_masks_zarr.get("masks").get("0"))

    # Get unique labels (skip 0 = background)
    labels = np.unique(nuclear_mask)
    labels = labels[labels > 0]

    # Get focus map and mean map
    focus_map = focus_maps.get("dapi_focus_map")
    mean_map = focus_maps.get("dapi_mean_map")

    if focus_map is None or mean_map is None:
        raise ValueError("focus_maps must contain 'dapi_focus_map' and 'dapi_mean_map'")

    # Aggregate pixel focus maps over nuclear masks using scipy.ndimage.mean
    logging.info("  Aggregating pixel focus maps over nuclear masks...")
    t0 = time.time()
    cell_focus = ndimage_mean(focus_map, nuclear_mask, labels)
    cell_intensity = ndimage_mean(mean_map, nuclear_mask, labels)
    logging.info(f"  [TIMING] ndimage.mean aggregation: {time.time() - t0:.1f}s")

    # Normalize: 99th percentile of mean intensity
    dapi_norm = np.percentile(cell_intensity, 99)
    if dapi_norm == 0:
        dapi_norm = 1.0
    ccfs_dapi = np.asarray(cell_focus) / dapi_norm

    # Get centroids and areas from regionprops (lightweight, no intensity needed)
    logging.info("  Computing nuclear centroids and areas...")
    t0 = time.time()
    from skimage.measure import regionprops

    props = regionprops(nuclear_mask)
    logging.info(f"  [TIMING] regionprops (centroids/areas): {time.time() - t0:.1f}s")

    # Build vectorized lookup from labels to CCFS/intensity values
    label_to_ccfs = dict(zip(labels, ccfs_dapi))
    label_to_intensity = dict(zip(labels, cell_intensity))

    # Build DataFrame from regionprops
    rows = []
    for p in props:
        lab = p.label
        cy, cx = p.centroid
        # Map nucleus centroid to cell ID
        iy = min(int(cy), cellseg_mask.shape[0] - 1)
        ix = min(int(cx), cellseg_mask.shape[1] - 1)
        cell_id = int(cellseg_mask[iy, ix])
        rows.append(
            {
                "label": lab,
                "centroid-0": cy,
                "centroid-1": cx,
                "area_nucleus": p.area,
                "CCFS_DAPI": label_to_ccfs.get(lab, np.nan),
                "mean_intensity": label_to_intensity.get(lab, np.nan),
                "CellID": cell_id,
            }
        )
    nucleus_props = pd.DataFrame(rows)
    del nuclear_mask, props

    # Cell-level areas from cellseg_mask (vectorized via bincount)
    logging.info("  Computing cell areas...")
    cell_labels_flat = cellseg_mask.ravel()
    cell_areas = np.bincount(cell_labels_flat)
    nucleus_props["area_cell"] = nucleus_props["CellID"].map(
        lambda cid: int(cell_areas[cid]) if cid < len(cell_areas) else 0
    )

    # Boundary and RNA per-cell intensities from pixel maps if available
    boundary_mean = focus_maps.get("boundary_mean_map")
    intrna_mean = focus_maps.get("intrna_mean_map")

    if boundary_mean is not None:
        logging.info("  Aggregating boundary intensity over cell masks...")
        cell_labels_unique = nucleus_props["CellID"].unique()
        cell_labels_unique = cell_labels_unique[cell_labels_unique > 0]
        boundary_per_cell = ndimage_mean(
            boundary_mean, cellseg_mask, cell_labels_unique
        )
        boundary_lookup = dict(zip(cell_labels_unique, boundary_per_cell))
        nucleus_props["mean_intensity_Boundary"] = nucleus_props["CellID"].map(
            boundary_lookup
        )
    else:
        nucleus_props["mean_intensity_Boundary"] = np.nan

    if intrna_mean is not None:
        logging.info("  Aggregating IntRNA intensity over cell masks...")
        cell_labels_unique = nucleus_props["CellID"].unique()
        cell_labels_unique = cell_labels_unique[cell_labels_unique > 0]
        intrna_per_cell = ndimage_mean(intrna_mean, cellseg_mask, cell_labels_unique)
        intrna_lookup = dict(zip(cell_labels_unique, intrna_per_cell))
        nucleus_props["mean_intensity_IntRNA"] = nucleus_props["CellID"].map(
            intrna_lookup
        )
    else:
        nucleus_props["mean_intensity_IntRNA"] = np.nan

    return nucleus_props


# ===== FALLBACK: REGIONPROPS-BASED CCFS (for legacy mode) =====


def calculate_ccfs_measurements(xoa_morphology_files, cellseg_mask, cell_masks_zarr):
    """
    Calculate CCFS measurements with improved variable naming and memory management.
    Loads data just before use and deletes it immediately after.
    """
    import numpy as np
    import pandas as pd

    # Load full resolution image channels only when needed
    fullres_channels = imread(
        xoa_morphology_files[0], is_ome=False, level=0, aszarr=False
    )

    # Check number of channels (could be 2D array for single channel, or 3D array for multi-channel)
    if len(fullres_channels.shape) == 2:
        # Single channel (2D array)
        dapi_image = fullres_channels
        boundary_image = None
        rna_image = None
    else:
        # Multi-channel (3D array: channels, height, width)
        dapi_image = fullres_channels[0]
        boundary_image = fullres_channels[1] if fullres_channels.shape[0] > 1 else None
        rna_image = fullres_channels[2] if fullres_channels.shape[0] > 2 else None
    del fullres_channels

    # Load nuclear segmentation mask
    nuclear_mask = np.array(cell_masks_zarr.get("masks").get("0"))

    # Per-nucleus measurements on DAPI image
    nucleus_props = pd.DataFrame(
        regionprops_table(dapi_image, nuclear_mask, position=True)
    )
    # Calculate CCFS_DAPI
    mean_intensity = nucleus_props["mean_intensity"]
    std_intensity = nucleus_props["standard_deviation_intensity"]
    dapi_norm = np.percentile(mean_intensity, 99)
    ccfs_dapi = (std_intensity * std_intensity) / (mean_intensity * dapi_norm)
    # Zero-intensity nuclei (mask over background) produce 0/0 = NaN.
    # No signal means no focus information: CCFS is 0.
    ccfs_dapi = ccfs_dapi.fillna(0.0)
    nucleus_props["CCFS_DAPI"] = ccfs_dapi
    del dapi_image

    # Get CellID for each nucleus by measuring mean intensity in cellseg_mask
    cellid_props = pd.DataFrame(
        regionprops_table(cellseg_mask, nuclear_mask, position=True)
    )
    nucleus_props["CellID"] = cellid_props["mean_intensity"].astype(int)
    del nuclear_mask, cellid_props

    # Measure boundary (red channel) intensity per cell
    if boundary_image is not None:
        boundary_props = pd.DataFrame(regionprops_table(boundary_image, cellseg_mask))
        boundary_props["CellID"] = boundary_props["label"].astype(int)
        del boundary_image
    else:
        # Create empty boundary_props if boundary channel not available
        boundary_props = pd.DataFrame(
            {"CellID": nucleus_props["CellID"], "mean_intensity": np.nan}
        )

    # Measure RNA (interior) intensity per cell
    if rna_image is not None:
        rna_props = pd.DataFrame(regionprops_table(rna_image, cellseg_mask))
        rna_props["CellID"] = rna_props["label"].astype(int)
        rna_props["mean_intensity_IntRNA"] = rna_props["mean_intensity"]
        del rna_image
    else:
        # Create empty rna_props if RNA channel not available
        rna_props = pd.DataFrame(
            {
                "CellID": nucleus_props["CellID"],
                "area": np.nan,
                "mean_intensity_IntRNA": np.nan,
            }
        )

    # Merge all measurements into a single DataFrame
    merged_data = pd.merge(
        nucleus_props,
        rna_props[["CellID", "area"]],
        how="inner",
        on="CellID",
        suffixes=("_nucleus", "_cell"),
    )
    merged_data = pd.merge(
        merged_data,
        boundary_props[["CellID", "mean_intensity"]],
        how="inner",
        on="CellID",
        suffixes=("_DAPI", "_Boundary"),
    )
    merged_data = pd.merge(
        merged_data,
        rna_props[["CellID", "mean_intensity_IntRNA"]],
        how="inner",
        on="CellID",
    )

    return merged_data


# ===== COMBINED MAIN =====


def _check_cell_data_exists(xenium_bundle_dir):
    """
    Check if cell data files exist in the Xenium bundle directory.

    Args:
        xenium_bundle_dir: Path to Xenium bundle directory

    Returns:
        bool: True if required cell data files exist
    """
    xenium_bundle_dir = Path(xenium_bundle_dir)
    required_files = [
        xenium_bundle_dir / "cells.parquet",
        xenium_bundle_dir / "cells.zarr.zip",
    ]
    # Check clustering path
    clusters_path = (
        xenium_bundle_dir
        / "analysis"
        / "clustering"
        / "gene_expression_kmeans_10_clusters"
        / "clusters.csv"
    )

    # Also check for analysis.tar.gz (test data)
    has_analysis = (
        clusters_path.exists() or (xenium_bundle_dir / "analysis.tar.gz").is_file()
    )

    for f in required_files:
        if not f.exists():
            logging.info(f"  Cell data check: {f.name} not found")
            return False

    if not has_analysis:
        logging.info("  Cell data check: clustering/UMAP data not found")
        return False

    return True


@click.command()
@click.option(
    "--xenium-bundle-dir", required=True, help="Path to Xenium bundle directory"
)
@click.option("--outdir", required=True, help="Output directory for results")
@click.option(
    "--stain-names",
    default=None,
    help="Semicolon-separated list of stain names. If not provided, uses defaults.",
)
@click.option(
    "--roi-size", default=35, type=int, show_default=True, help="Tile size in pixels"
)
@click.option(
    "--max-scatter-points",
    default=10000,
    type=int,
    show_default=True,
    help="Maximum number of points to plot in scatter figures. Set to 0 to plot all points.",
)
@click.option(
    "--legacy-focus",
    is_flag=True,
    default=False,
    help="Use legacy per-tile loop focus scoring instead of the default convolution-based GPU-accelerated method.",
)
@click.option(
    "--sample-id",
    default=None,
    help="Sample identifier for logging and metrics output.",
)
@click.option(
    "--no-snr",
    is_flag=True,
    default=False,
    help="Disable SNR metrics (image Otsu / quartiles, transcripts, slide matrix, neg spatial).",
)
@click.option(
    "--snr-no-roi-tx-table",
    is_flag=True,
    default=False,
    help="Do not write SNR_roi_tx.parquet (or .csv.gz) alongside roi_qc_metrics.",
)
@click.option(
    "--snr-otsu-max-rois",
    type=int,
    default=None,
    help="Cap tiles for per-tile Otsu image SNR (default: all tiles).",
)
@click.option(
    "--snr-with-moran",
    is_flag=True,
    default=False,
    help="SNR only: enable Moran's I for neg spatial (needs PySAL/esda; default off).",
)
@click.option(
    "--save-dapi-maps-tiff",
    "save_dapi_maps_tiff",
    is_flag=True,
    default=False,
    help=(
        "Write full-resolution per-pixel maps as tiled float32 TIFF "
        "(dapi_focus/mean/lap_var and boundary/intrna focus/mean when present). "
        "Default off — QC and SNR use in-memory arrays only; enable for archival "
        "or external tools (large files)."
    ),
)
@click.option(
    "--roi-thresholds-yaml",
    default=None,
    type=click.Path(exists=True),
    help="Path to tile image QC thresholds YAML. Overrides hardcoded defaults.",
)
@click.option(
    "--lap-sigma",
    default=1.0,
    type=float,
    show_default=True,
    help="Gaussian sigma for Laplacian of Gaussian (LoG) pre-smoothing.",
)
def main(
    xenium_bundle_dir,
    outdir,
    stain_names,
    roi_size,
    max_scatter_points,
    legacy_focus,
    sample_id,
    no_snr,
    snr_no_roi_tx_table,
    snr_otsu_max_rois,
    snr_with_moran,
    save_dapi_maps_tiff,
    roi_thresholds_yaml,
    lap_sigma,
):
    """
    Combined Xenium Image QC pipeline.

    Performs pixel-level focus analysis (GPU-accelerated), tile-level analysis,
    and optionally cell-level analysis when cell data is available.

    Produces: tile figures, cell figures, image_qc_metrics.json, image_qc_metrics.csv,
    13 Quarto-required PNGs, and versions.yml.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    global _MAX_SCATTER_POINTS
    if max_scatter_points > 0:
        _MAX_SCATTER_POINTS = max_scatter_points

    t_total_start = time.time()
    _timings: dict = {}  # Collect per-step wall-clock seconds for profiling
    logging.info("=" * 60)
    logging.info("Starting Combined Xenium Image QC Pipeline")
    if sample_id:
        logging.info("Sample ID: %s", sample_id)
    logging.info("=" * 60)
    logging.info(f"Input directory: {xenium_bundle_dir}")
    logging.info(f"Output directory: {outdir}")

    # Load QC thresholds from YAML (or use hardcoded defaults)
    qc_thresholds = _load_qc_thresholds(roi_thresholds_yaml)
    if roi_thresholds_yaml:
        logging.info(f"Loaded QC thresholds from: {roi_thresholds_yaml}")
    else:
        logging.info("Using built-in default QC thresholds")
    # Resolve lap_sigma from YAML (CLI value is the fallback)
    _yaml_sigma = qc_thresholds.get("lap_sigma", lap_sigma)
    try:
        _lap_sigma = float(_yaml_sigma)
    except (TypeError, ValueError):
        logging.warning(
            "Invalid lap_sigma value %r in YAML, falling back to CLI=%s",
            _yaml_sigma,
            lap_sigma,
        )
        _lap_sigma = float(lap_sigma)
    logging.info(
        f"LoG sigma: {_lap_sigma} (CLI={lap_sigma}, YAML={qc_thresholds.get('lap_sigma', 'not set')})"
    )

    # Resolve top-level operational thresholds from YAML (with hardcoded fallbacks)
    _roi_intensity_threshold = float(
        qc_thresholds.get("roi_intensity_threshold", ROI_INTENSITY_THRESHOLD)
    )
    _min_tissue_cov = float(
        qc_thresholds.get(
            "min_tissue_coverage_for_intensity_qc",
            ROI_MIN_TISSUE_COVERAGE_FOR_INTENSITY_QC,
        )
    )
    _roi_focus_pct = float(
        qc_thresholds.get("roi_focus_score_percentile", ROI_FOCUS_SCORE_PERCENTILE)
    )
    _blur_prob_thresh = float(qc_thresholds.get("blur_prob_threshold", 0.5))
    _focus_cfg = qc_thresholds.get("focus") or {}
    _ccfs_low_texture_threshold = float(
        _focus_cfg.get("ccfs_low_texture_threshold", DEFAULT_CCFS_LOW_TEXTURE_THRESHOLD)
    )
    logging.info(
        f"  roi_intensity_threshold={_roi_intensity_threshold}, "
        f"min_tissue_cov={_min_tissue_cov}, "
        f"roi_focus_pct={_roi_focus_pct}, "
        f"blur_prob_thresh={_blur_prob_thresh}, "
        f"ccfs_low_texture_threshold={_ccfs_low_texture_threshold}"
    )

    # ===== PHASE 1: LOAD =====
    logging.info("\n--- PHASE 1: LOAD ---")

    # Handle stain_names (default if None)
    if stain_names is None:
        stain_names_list = [
            "DAPI",
            "Boundary (ATP1A1/E-Cadherin/CD45)",
            "Interior - RNA (18S)",
            "Protein (alphaSMA/Vimentin)",
        ]
        logging.info(f"Using default stain names: {stain_names_list}")
    else:
        logging.info(f"Stain names: {stain_names}")
        if ";" in stain_names:
            stain_names_list = stain_names.split(";")
        else:
            stain_names_list = [stain_names]

    # Load downsampled images for sample QC
    morphology_focus_dir = Path(xenium_bundle_dir) / "morphology_focus"
    if (morphology_focus_dir / "morphology_focus_0000.ome.tif").exists():
        xoa_morphology_files = [
            morphology_focus_dir / "morphology_focus_0000.ome.tif",
            morphology_focus_dir / "morphology_focus_0001.ome.tif",
            morphology_focus_dir / "morphology_focus_0002.ome.tif",
            morphology_focus_dir / "morphology_focus_0003.ome.tif",
        ]
    else:
        xoa_morphology_files = sorted(
            list(morphology_focus_dir.glob("ch000*.ome.tif")),
            key=lambda x: x.stem.split("_")[0],
        )

    if len(xoa_morphology_files) < 1 or not xoa_morphology_files[0].exists():
        raise ValueError(
            f"Expected at least 1 morphology file (DAPI), got "
            f"{len(xoa_morphology_files)}"
        )

    # Load and prepare data (creates output directories)
    t0 = time.time()
    data = load_and_prepare_data(xenium_bundle_dir, outdir)
    _timings["load_and_prepare_data"] = time.time() - t0

    # Load morphology images (level 3, downsampled)
    logging.info("Loading morphology images...")
    t0 = time.time()
    small0, small1, small2 = load_morphology_images(xoa_morphology_files)
    logging.info(
        f"Loaded morphology images: {small0.shape}, "
        f"{small1.shape if small1 is not None else 'None'}, "
        f"{small2.shape if small2 is not None else 'None'}"
    )
    _timings["load_morphology_images"] = time.time() - t0
    logging.info(
        f"[TIMING] Loading morphology images: {_timings['load_morphology_images']:.1f}s"
    )

    # Check if cell data exists
    has_cell_data = _check_cell_data_exists(xenium_bundle_dir)
    if has_cell_data:
        logging.info("Cell data detected - will perform cell-level analysis")
    else:
        logging.info("No cell data detected - will skip cell-level analysis")

    # ===== PHASE 2: PIXEL-LEVEL FOCUS MAPS =====
    logging.info("\n--- PHASE 2: PIXEL-LEVEL FOCUS MAPS ---")

    # Generate tissue masks and distance maps (cell-independent)
    logging.info("Generating tissue masks and distance maps...")
    t0 = time.time()
    whole_sample, holes, dense_intensity_regions, distance_map, distance_map2 = (
        generate_tissue_mask(xoa_morphology_files, small0, small1, small2)
    )
    _timings["generate_tissue_mask"] = time.time() - t0
    logging.info("Generated tissue masks and distance maps")
    logging.info(
        f"[TIMING] Tissue mask generation: {_timings['generate_tissue_mask']:.1f}s"
    )

    # Validate ROI size
    if roi_size <= 0:
        roi_size = 35
        logging.info(f"Using default tile size: {roi_size}px")
    else:
        logging.info(f"Using tile size: {roi_size}px")

    # Auto-detect available GPUs
    available_gpus = detect_gpu_ids()
    if available_gpus:
        logging.info(f"Detected {len(available_gpus)} GPU(s): {available_gpus}")
    else:
        logging.info("No GPUs detected, using CPU backend")

    # Calculate cell-independent grid ROI focus scores
    logging.info("Calculating cell-independent grid tile focus scores...")
    t0 = time.time()
    if legacy_focus:
        logging.info("  Using LEGACY per-tile loop method (--legacy-focus)")
        df_grid_roi = calculate_roi_focusscore_without_laplace(
            xoa_morphology_files,
            roi_size=roi_size,
            stride=None,
            tissue_filter=True,
            min_tissue_coverage=0.0,
        )
        focus_maps = None
    else:
        logging.info("  Using convolution-based GPU-accelerated method")
        df_grid_roi, focus_maps = calculate_roi_focusscore(
            xoa_morphology_files,
            roi_size=roi_size,
            stride=None,
            tissue_filter=True,
            min_tissue_coverage=0.0,
            gpu_ids=available_gpus if available_gpus else None,
            return_pixel_maps=True,
            lap_sigma=_lap_sigma,
        )
    _timings["calculate_roi_focusscore"] = time.time() - t0
    logging.info(f"Calculated grid tile focus scores for {len(df_grid_roi):,} tiles")
    logging.info(
        f"[TIMING] calculate_roi_focusscore(): {_timings['calculate_roi_focusscore']:.1f}s"
    )

    # ===== PHASE 3: ROI-LEVEL ANALYSIS =====
    logging.info("\n--- PHASE 3: TILE-LEVEL ANALYSIS ---")

    # Calculate ROI blur threshold
    logging.info("Calculating tile blur threshold...")
    t0 = time.time()
    roi_threshold = calculate_roi_blur_threshold(
        df_grid_roi,
        intensity_threshold=_roi_intensity_threshold,
        focus_percentile=_roi_focus_pct,
    )
    _timings["calculate_roi_blur_threshold"] = time.time() - t0
    logging.info(f"  Calculated threshold: {roi_threshold:.2f} (raw score units)")
    logging.info(f"  Intensity threshold: {_roi_intensity_threshold}")

    # Fit 1D GMM model
    logging.info("Fitting 1D GMM model for tile focus scores...")
    t0 = time.time()
    gmm, blur_component_idx = fit_focus_gmm(
        df_grid_roi,
        intensity_threshold=_roi_intensity_threshold,
        focus_col_name="dapi_focus_score",
    )

    logging.info("Classifying tiles (1D GMM)...")
    df_grid_roi = classify_roi_blur(
        df_grid_roi,
        gmm=gmm,
        blur_component_idx=blur_component_idx,
        blur_prob_threshold=_blur_prob_thresh,
        intensity_threshold=_roi_intensity_threshold,
        focus_col_name="dapi_focus_score",
    )
    n_blurred = int(df_grid_roi["is_blurred_gmm"].sum())
    n_in_focus = int((~df_grid_roi["is_blurred_gmm"]).sum())
    _timings["fit_gmm_and_classify_1d"] = time.time() - t0
    logging.info(f"  Blurred: {n_blurred:,}, In-focus: {n_in_focus:,} (1D GMM)")
    logging.info(
        f"[TIMING] GMM fit + classify (1D): {_timings['fit_gmm_and_classify_1d']:.1f}s"
    )

    # Fit 2D GMM if Laplacian variance available
    gmm_2d = None
    blur_component_idx_2d = None
    t0 = time.time()
    if "dapi_lap_var" in df_grid_roi.columns:
        logging.info("\nFitting 2D GMM model (focus_score + Laplacian variance)...")
        try:
            gmm_2d, blur_component_idx_2d = fit_focus_gmm_2d(
                df_grid_roi,
                intensity_threshold=_roi_intensity_threshold,
                focus_col_name="dapi_focus_score",
            )
            logging.info("Classifying tiles (2D GMM)...")
            df_grid_roi = classify_roi_blur_2d(
                df_grid_roi,
                gmm=gmm_2d,
                blur_component_idx=blur_component_idx_2d,
                blur_prob_threshold=_blur_prob_thresh,
                intensity_threshold=_roi_intensity_threshold,
                focus_col_name="dapi_focus_score",
            )
            n_blurred_2d = int(df_grid_roi["is_blurred_gmm_2d"].sum())
            n_in_focus_2d = int((~df_grid_roi["is_blurred_gmm_2d"]).sum())
            logging.info(
                f"  Blurred: {n_blurred_2d:,}, In-focus: {n_in_focus_2d:,} (2D GMM)"
            )
        except Exception as e:
            logging.warning(f"  Warning: 2D GMM failed: {e}")
            logging.info("  Continuing with 1D GMM results only")
            gmm_2d = None
            blur_component_idx_2d = None
    else:
        logging.info("\nSkipping 2D GMM: dapi_lap_var column not found")
    _timings["fit_gmm_2d_and_classify"] = time.time() - t0
    logging.info(
        f"[TIMING] GMM 2D fit + classify: {_timings['fit_gmm_2d_and_classify']:.1f}s"
    )

    # Save grid ROI data to CSV
    t0 = time.time()
    grid_roi_csv = data["outdir"] / "grid_roi_focus_scores.csv"
    df_grid_roi.to_csv(grid_roi_csv, index=False)
    _timings["save_grid_roi_csv"] = time.time() - t0
    logging.info(f"Saved grid tile data to {grid_roi_csv}")

    # Optional: write full-res per-pixel maps (very large); not needed for QC/SNR in-process
    if focus_maps is not None and save_dapi_maps_tiff:
        logging.info("Saving pixel-level focus maps as TIFF (--save-dapi-maps-tiff)...")
        t0 = time.time()
        save_pixel_focus_maps(focus_maps, data["outdir"])
        _timings["save_pixel_focus_maps"] = time.time() - t0
        logging.info(
            f"[TIMING] save_pixel_focus_maps(): {_timings['save_pixel_focus_maps']:.1f}s"
        )
    elif focus_maps is not None:
        logging.info(
            "Skipping pixel-map TIFF export (default). "
            "Pass --save-dapi-maps-tiff to write dapi_*/boundary_*/intrna_* map TIFFs."
        )

    # Free Laplacian map — no longer needed after TIFF export / GMM fitting.
    # Remaining consumers (SNR, figures, CCFS) only use focus/mean maps.
    if focus_maps is not None:
        for _k in ("dapi_lap_var_map",):
            focus_maps.pop(_k, None)

    # Save threshold configuration
    total_rois_gmm = int(len(df_grid_roi))
    n_blurred_gmm = int(df_grid_roi["is_blurred_gmm"].sum())
    pct_blurred_gmm = (
        (n_blurred_gmm / total_rois_gmm * 100.0) if total_rois_gmm > 0 else 0.0
    )

    threshold_config = {
        "roi_focus_score_threshold": float(roi_threshold),
        "roi_intensity_threshold": float(_roi_intensity_threshold),
        "roi_focus_score_percentile": float(_roi_focus_pct),
        "threshold_method": "Option B: Percentile from tissue tiles (intensity >= threshold) + fixed intensity threshold",
        "gmm_1d": {
            "n_components": int(gmm.n_components),
            "blur_component_index": int(blur_component_idx),
            "component_means_log1p_focus": [float(m) for m in gmm.means_.flatten()],
            "component_weights": [float(w) for w in gmm.weights_.flatten()],
            "blur_prob_threshold": float(_blur_prob_thresh),
            "fraction_rois_blurred_gmm": float(pct_blurred_gmm),
            "features": ["log1p(dapi_focus_score)"],
        },
        "units": {
            "roi_focus_score_threshold_percentile": "raw score units",
            "roi_intensity_threshold": "raw pixel intensity (16-bit, 0-65535)",
            "component_means_log1p_focus": "log1p(raw focus score)",
        },
    }

    if "is_blurred_gmm_2d" in df_grid_roi.columns and gmm_2d is not None:
        n_blurred_gmm_2d = int(df_grid_roi["is_blurred_gmm_2d"].sum())
        pct_blurred_gmm_2d = (
            (n_blurred_gmm_2d / total_rois_gmm * 100.0) if total_rois_gmm > 0 else 0.0
        )
        threshold_config["gmm_2d"] = {
            "n_components": int(gmm_2d.n_components),
            "blur_component_index": int(blur_component_idx_2d),
            "component_means": [[float(m[0]), float(m[1])] for m in gmm_2d.means_],
            "component_weights": [float(w) for w in gmm_2d.weights_.flatten()],
            "blur_prob_threshold": float(_blur_prob_thresh),
            "fraction_rois_blurred_gmm": float(pct_blurred_gmm_2d),
            "features": ["log1p(dapi_focus_score)", "log1p(dapi_lap_var)"],
        }
        threshold_config["units"]["component_means_2d"] = (
            "[log1p(focus_score), log1p(lap_var)]"
        )

    threshold_json = data["outdir"] / "roi_blur_threshold.json"
    with open(threshold_json, "w") as f:
        json.dump(threshold_config, f, indent=2)
    logging.info(f"Saved tile blur threshold configuration to {threshold_json}")

    # Save ROI count
    roi_count_file = data["outdir"] / "roi_count.txt"
    with open(roi_count_file, "w") as f:
        f.write(str(len(df_grid_roi)))

    # Calculate ROI intensities
    logging.info("Calculating tile intensities...")
    t0 = time.time()
    df_roi_intensities = calculate_roi_intensities(xoa_morphology_files, df_grid_roi)
    _timings["calculate_roi_intensities"] = time.time() - t0
    logging.info(
        f"[TIMING] Tile intensity calculation: {_timings['calculate_roi_intensities']:.1f}s"
    )

    snr_summary = None
    if not no_snr:
        logging.info("Computing SNR metrics (roi_qc / snr_metrics)...")
        t_snr = time.time()
        pix_um = snr_metrics.read_xenium_pixel_size_um(Path(xenium_bundle_dir))
        if pix_um is None:
            logging.info(
                "  experiment.xenium missing or invalid pixel_size — transcript SNR may skip"
            )
        try:
            df_roi_intensities, snr_summary = snr_metrics.compute_snr_summary(
                df_roi_intensities,
                bundle_dir=Path(xenium_bundle_dir),
                outdir=data["outdir"],
                focus_maps=focus_maps,
                pixel_size_um=pix_um,
                otsu_max_rois=snr_otsu_max_rois,
                save_roi_tx_table=not snr_no_roi_tx_table,
                snr_include_moran=snr_with_moran,
                # Same stride as grid in calculate_roi_focusscore (stride=None → roi_size).
                roi_grid_stride=(roi_size, roi_size),
                snr_thresholds=qc_thresholds.get("snr") or {},
            )
        except Exception as e:
            logging.warning("SNR metrics failed (continuing image QC): %s", e)
            snr_summary = {
                "status": "error",
                "error": str(e),
                "components": {},
                "verdict": {"overall_snr_verdict": "NOT_COMPUTED"},
            }
        _timings["snr_metrics"] = time.time() - t_snr
        logging.info(f"[TIMING] SNR metrics: {_timings['snr_metrics']:.1f}s")
    else:
        logging.info("SNR metrics disabled (--no-snr)")
    df_grid_roi = df_roi_intensities

    # Assess raw intensity quality (thresholds from YAML or defaults)
    _ch_cfg = qc_thresholds.get("channels") or {}
    _ic_dapi = (_ch_cfg.get("DAPI") or {}).get(
        "intensity_critical", _INTENSITY_CRITICAL_DEFAULTS["dapi"]
    )
    _ic_boundary = (_ch_cfg.get("boundary") or {}).get(
        "intensity_critical", _INTENSITY_CRITICAL_DEFAULTS["boundary"]
    )
    _ic_intrna = (_ch_cfg.get("intRNA") or {}).get(
        "intensity_critical", _INTENSITY_CRITICAL_DEFAULTS["intrna"]
    )
    logging.info("Assessing raw intensity quality...")
    t0 = time.time()
    intensity_stats = assess_raw_intensity_quality(
        df_roi_intensities,
        dapi_threshold_critical=_ic_dapi,
        boundary_threshold_critical=_ic_boundary,
        intrna_threshold_critical=_ic_intrna,
        min_tissue_coverage=_min_tissue_cov,
        channel_pct_thresholds=_ch_cfg,
    )
    _timings["assess_raw_intensity_quality"] = time.time() - t0
    logging.info(
        f"[TIMING] assess_raw_intensity_quality(): {_timings['assess_raw_intensity_quality']:.1f}s"
    )

    # Save intensity statistics
    intensity_json = data["outdir"] / "intensity_assessment.json"
    with open(intensity_json, "w") as f:
        json.dump(intensity_stats, f, indent=2)

    # Generate ROI figures (cell-independent)
    logging.info("Generating tile-level figures...")
    t0 = time.time()
    generate_roi_figures(
        data,
        small0,
        small1,
        small2,
        distance_map,
        distance_map2,
        whole_sample,
        holes,
        dense_intensity_regions,
        df_grid_roi,
        df_roi_intensities,
        intensity_stats,
        xoa_morphology_files,
        focus_maps=focus_maps,
        snr_thresholds=(qc_thresholds.get("snr") or {}).get("roi_tx") or {},
        blur_prob_threshold=_blur_prob_thresh,
    )
    _timings["generate_roi_figures"] = time.time() - t0
    logging.info(
        f"[TIMING] generate_roi_figures(): {_timings['generate_roi_figures']:.1f}s"
    )

    # Free focus-score maps no longer needed.  CCFS only requires
    # dapi_focus_map, dapi_mean_map, boundary_mean_map, intrna_mean_map.
    if focus_maps is not None:
        for _k in ("boundary_focus_map", "intrna_focus_map"):
            focus_maps.pop(_k, None)

    # Save ROI QC metrics
    logging.info("Saving tile QC metrics...")
    t0 = time.time()
    save_roi_qc_metrics(
        df_grid_roi,
        intensity_stats,
        data["outdir"],
        roi_size=roi_size,
        snr_summary=snr_summary,
        distance_map=distance_map,
        distance_map2=distance_map2,
        edge_distance_threshold=-25.0,
        hole_distance_threshold=-25.0,
        min_tissue_coverage_for_qc=_min_tissue_cov,
        qc_thresholds=qc_thresholds,
        lap_sigma=_lap_sigma,
    )
    _timings["save_roi_qc_metrics"] = time.time() - t0
    logging.info(
        f"[TIMING] save_roi_qc_metrics(): {_timings['save_roi_qc_metrics']:.1f}s"
    )

    # ===== PHASE 4: CELL-LEVEL ANALYSIS (conditional on cell data) =====
    if has_cell_data:
        logging.info("\n--- PHASE 4: CELL-LEVEL ANALYSIS ---")

        # Prepare cell-centred output directory
        figures_dir_cell = data["outdir"] / "figures"
        figures_dir_cell.mkdir(parents=True, exist_ok=True)
        figures_source_dir_cell = figures_dir_cell / "figures_source"
        figures_source_dir_cell.mkdir(parents=True, exist_ok=True)

        # Prepare data dict for cell-level functions
        cell_data = dict(data)
        cell_data["figures_dir"] = figures_dir_cell

        # Unpack analysis.tar.gz if needed (test data)
        xenium_bundle_dir_path = Path(xenium_bundle_dir)
        if (xenium_bundle_dir_path / "analysis.tar.gz").is_file():
            import shutil

            shutil.unpack_archive(
                xenium_bundle_dir_path / "analysis.tar.gz", extract_dir="."
            )
            cell_data["clusters_csv_path"] = (
                Path("analysis")
                / "clustering"
                / "gene_expression_kmeans_10_clusters"
                / "clusters.csv"
            )
            cell_data["umap_path"] = (
                Path("analysis")
                / "umap"
                / "gene_expression_2_components"
                / "projection.csv"
            )
        else:
            cell_data["clusters_csv_path"] = (
                xenium_bundle_dir_path
                / "analysis"
                / "clustering"
                / "gene_expression_kmeans_10_clusters"
                / "clusters.csv"
            )
            cell_data["umap_path"] = (
                xenium_bundle_dir_path
                / "analysis"
                / "umap"
                / "gene_expression_2_components"
                / "projection.csv"
            )
        cell_data["cells_parquet_path"] = xenium_bundle_dir_path / "cells.parquet"
        cell_data["cell_masks_path"] = xenium_bundle_dir_path / "cells.zarr.zip"

        # Load spatial data
        logging.info("Loading spatial data and masks...")
        t0 = time.time()
        try:
            df_spatial, cellseg_mask, cell_masks_zarr = load_spatial_data(cell_data)
            logging.info(f"Loaded {len(df_spatial):,} cells")
            _timings["load_spatial_data"] = time.time() - t0
            logging.info(
                f"[TIMING] Loading spatial data: {_timings['load_spatial_data']:.1f}s"
            )
        except (FileNotFoundError, KeyError, AttributeError) as e:
            logging.warning(f"Warning: Could not load spatial data: {e}")
            logging.info("Skipping cell-level analysis")
            has_cell_data = False

    if has_cell_data:
        # Map distances to cells
        logging.info("Mapping distances to cells...")
        t0 = time.time()
        df_spatial = map_distances_to_cells(
            df_spatial, distance_map, distance_map2, dense_intensity_regions
        )
        _timings["map_distances_to_cells"] = time.time() - t0
        logging.info(
            f"[TIMING] map_distances_to_cells(): {_timings['map_distances_to_cells']:.1f}s"
        )

        # Calculate CCFS measurements
        logging.info("Calculating CCFS measurements...")
        t0 = time.time()
        if (
            focus_maps is not None
            and "dapi_focus_map" in focus_maps
            and "dapi_mean_map" in focus_maps
        ):
            # NEW: Use pixel-map aggregation
            logging.info("  Using pixel-map aggregation (scipy.ndimage.mean)")
            myData = calculate_ccfs_from_focus_maps(
                focus_maps, cell_masks_zarr, cellseg_mask, xoa_morphology_files
            )
        else:
            # Fallback: use regionprops-based method
            logging.info("  Using regionprops-based method (legacy)")
            myData = calculate_ccfs_measurements(
                xoa_morphology_files, cellseg_mask, cell_masks_zarr
            )
        _timings["ccfs_calculation"] = time.time() - t0
        logging.info(f"Calculated CCFS for {len(myData):,} cells")
        logging.info(f"[TIMING] CCFS calculation: {_timings['ccfs_calculation']:.1f}s")

        # All pixel-level focus maps consumed — free remaining memory.
        del focus_maps
        focus_maps = None

        # Add boolean columns to myData for thresholding (needed by generate_all_figures)
        ccfs_threshold = _ccfs_low_texture_threshold
        myData["is_low_nuclear_texture"] = myData["CCFS_DAPI"] <= ccfs_threshold
        myData["is_high_nuclear_texture"] = myData["CCFS_DAPI"] > ccfs_threshold

        # Map ROI focus scores to cells
        logging.info("Mapping tile focus scores to cells...")
        t0 = time.time()
        roi_mapped = map_grid_roi_to_cells(df_grid_roi, df_spatial, overlapping=False)
        logging.info(
            f"Mapped tile focus scores to {len(roi_mapped[roi_mapped['DAPI_RFSnorm_roi'].notna()]):,} cells"
        )
        _timings["map_grid_roi_to_cells"] = time.time() - t0
        logging.info(
            f"[TIMING] map_grid_roi_to_cells(): {_timings['map_grid_roi_to_cells']:.1f}s"
        )

        # Load ROI blur threshold
        roi_threshold_cell, roi_intensity_threshold = load_roi_blur_threshold(
            data["outdir"]
        )
        if roi_threshold_cell is None or roi_intensity_threshold is None:
            roi_threshold_cell = roi_threshold
            roi_intensity_threshold = _roi_intensity_threshold

        # Create merged dataset (superset version with ROI data)
        logging.info("Creating final merged dataset...")
        t0 = time.time()
        new_df, calculated_roi_threshold = create_final_merged_data(
            df_spatial,
            myData,
            roi_data=roi_mapped,
            ccfs_threshold=_ccfs_low_texture_threshold,
            roi_threshold=roi_threshold_cell,
            roi_intensity_threshold=roi_intensity_threshold,
        )
        _timings["create_final_merged_data"] = time.time() - t0
        logging.info(f"Created merged dataset with {len(new_df):,} cells")
        logging.info(
            f"[TIMING] create_final_merged_data(): {_timings['create_final_merged_data']:.1f}s"
        )

        # Alias dense_intensity_regions as artefacts for generate_all_figures compatibility
        artefacts = dense_intensity_regions

        # Ensure In-Area-with-Artefact column exists (needed by generate_all_figures)
        if "In-Area-with-Artefact" not in df_spatial.columns:
            # Map dense_intensity_regions to the artefact column name
            if "Dense-Intensity-Region-ID" in df_spatial.columns:
                df_spatial["In-Area-with-Artefact"] = df_spatial[
                    "Dense-Intensity-Region-ID"
                ]
            else:
                df_spatial["In-Area-with-Artefact"] = 0

        if "In-Area-with-Artefact" not in new_df.columns:
            if "Dense-Intensity-Region-ID" in new_df.columns:
                new_df["In-Area-with-Artefact"] = new_df["Dense-Intensity-Region-ID"]
            else:
                new_df["In-Area-with-Artefact"] = 0

        # Ensure has_artifacts column exists
        if "has_artifacts" not in new_df.columns:
            new_df["has_artifacts"] = new_df.get(
                "has_dense_intensity_regions",
                new_df.get("In-Area-with-Artefact", 0) > 0,
            )

        # Generate 13 Quarto-required figures
        logging.info("Generating Quarto-required figures (13 PNGs)...")
        t0 = time.time()
        generate_all_figures(
            cell_data,
            df_spatial,
            new_df,
            myData,
            small0,
            small1,
            small2,
            distance_map,
            distance_map2,
            whole_sample,
            holes,
            artefacts,
            ccfs_low_texture_threshold=_ccfs_low_texture_threshold,
        )
        _timings["generate_all_figures"] = time.time() - t0
        logging.info(
            f"[TIMING] generate_all_figures(): {_timings['generate_all_figures']:.1f}s"
        )

        # Generate cell-centred comparison figures (ROI vs CCFS)
        figures_cell_centred_dir = data["outdir"] / "figures_cell_centred"
        figures_cell_centred_dir.mkdir(parents=True, exist_ok=True)
        figures_cell_centred_source = figures_cell_centred_dir / "figures_source"
        figures_cell_centred_source.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        generate_cell_figures(
            cell_data,
            new_df,
            myData,
            figures_cell_centred_dir,
            figures_cell_centred_source,
            roi_threshold=calculated_roi_threshold,
            roi_intensity_threshold=roi_intensity_threshold,
            ccfs_low_texture_threshold=_ccfs_low_texture_threshold,
        )
        _timings["generate_cell_figures"] = time.time() - t0
        logging.info(
            f"[TIMING] generate_cell_figures(): {_timings['generate_cell_figures']:.1f}s"
        )

        # Save cell QC metrics (superset version with ROI metrics)
        logging.info("Saving cell QC metrics...")
        t0 = time.time()
        save_cell_qc_metrics(new_df, data["outdir"], roi_size=roi_size)
        _timings["save_cell_qc_metrics"] = time.time() - t0
        logging.info(
            f"[TIMING] save_cell_qc_metrics(): {_timings['save_cell_qc_metrics']:.1f}s"
        )

        # Save simple image_qc_metrics.json for Quarto compatibility
        logging.info("Saving Quarto-compatible image_qc_metrics.json...")
        n_total = len(new_df)
        n_ccfs_low_texture = int(new_df["is_low_nuclear_texture"].sum())
        qc_metrics = {
            "total_cells": n_total,
            "cells_with_transcripts": len(new_df[new_df["transcript_counts"] > 0]),
            "mean_transcript_count": float(new_df["transcript_counts"].mean()),
            "median_transcript_count": float(new_df["transcript_counts"].median()),
            "mean_ccfs_dapi": float(new_df["CCFS_DAPI"].mean()),
            "median_ccfs_dapi": float(new_df["CCFS_DAPI"].median()),
            "ccfs_low_texture_threshold": float(_ccfs_low_texture_threshold),
            "cells_high_nuclear_texture": len(
                new_df[new_df["is_high_nuclear_texture"]]
            ),
            "cells_low_nuclear_texture": n_ccfs_low_texture,
            "pct_low_nuclear_texture": round(100.0 * n_ccfs_low_texture / n_total, 4)
            if n_total > 0
            else 0.0,
            "cells_near_edge": len(new_df[new_df["is_near_edge"]]),
            "cells_near_holes": len(new_df[new_df["is_near_hole"]]),
            "cells_with_artifacts": int(new_df["has_artifacts"].sum()),
            "clusters_present": sorted(new_df["Cluster_kmeans10"].unique().tolist()),
            "segmentation_methods": sorted(
                new_df["segmentation_method"].unique().tolist()
            ),
        }
        # GMM-ROI blur metrics (if cell-to-ROI mapping was performed)
        if "is_blurred_gmm_2d_roi" in new_df.columns:
            n_gmm = int(new_df["is_blurred_gmm_2d_roi"].sum())
            qc_metrics["cells_blurred_gmm_2d_roi"] = n_gmm
            qc_metrics["pct_blurred_gmm_2d_roi"] = (
                round(100.0 * n_gmm / n_total, 4) if n_total > 0 else 0.0
            )
            agreement = int(
                (
                    new_df["is_low_nuclear_texture"] == new_df["is_blurred_gmm_2d_roi"]
                ).sum()
            )
            qc_metrics["ccfs_gmm_agreement_pct"] = (
                round(100.0 * agreement / n_total, 4) if n_total > 0 else 0.0
            )
        # Per-cluster blur stats for cluster-outlier detection
        if "Cluster_kmeans10" in new_df.columns:
            blur_col = (
                "is_blurred_gmm_2d_roi"
                if "is_blurred_gmm_2d_roi" in new_df.columns
                else "is_low_nuclear_texture"
            )
            cluster_blur = {}
            for cl, grp in new_df.groupby("Cluster_kmeans10"):
                n_cl = len(grp)
                n_blur_cl = int(grp[blur_col].sum())
                pct_blur_cl = round(100.0 * n_blur_cl / n_cl, 2) if n_cl > 0 else 0.0
                cluster_blur[int(cl)] = {
                    "n_cells": n_cl,
                    "n_blurred": n_blur_cl,
                    "pct_blurred": pct_blur_cl,
                    "median_ccfs_dapi": round(float(grp["CCFS_DAPI"].median()), 4)
                    if not pd.isna(grp["CCFS_DAPI"].median())
                    else 0.0,
                }
            qc_metrics["cluster_blur"] = cluster_blur
            qc_metrics["cluster_blur_method"] = blur_col
            # Identify outlier clusters (>2× sample-wide blur rate)
            sample_pct = qc_metrics.get(
                "pct_blurred_gmm_2d_roi", qc_metrics.get("pct_low_nuclear_texture", 0.0)
            )
            outlier_threshold = max(
                sample_pct * 2.0, 10.0
            )  # at least 10% to avoid noise
            outlier_clusters = {
                k: v
                for k, v in cluster_blur.items()
                if v["pct_blurred"] > outlier_threshold and v["n_cells"] >= 50
            }
            if outlier_clusters:
                qc_metrics["cluster_blur_outliers"] = outlier_clusters

        with open(data["outdir"] / "image_qc_metrics.json", "w") as f:
            json.dump(qc_metrics, f, indent=2)

        # Save image_qc_metrics.csv
        cell_id_col = new_df.pop("cell_id")
        new_df.insert(0, "cell_id", cell_id_col)
        new_df.to_csv(data["outdir"] / "image_qc_metrics.csv", index=False)

        # Save dense intensity region summary
        if "Dense-Intensity-Region-ID" in new_df.columns:
            region_summary = (
                new_df[new_df["Dense-Intensity-Region-ID"] > 0]
                .groupby("Dense-Intensity-Region-ID")
                .agg({"cell_id": "count", "x": "mean", "y": "mean"})
                .reset_index()
            )
            region_summary.columns = [
                "Dense-Intensity-Region-ID",
                "cell_count",
                "mean_x",
                "mean_y",
            ]
            region_summary = region_summary.sort_values("Dense-Intensity-Region-ID")
            region_summary["annotation"] = ""
            region_summary.to_csv(
                data["outdir"] / "dense_intensity_regions_summary.csv", index=False
            )
    else:
        logging.info("\n--- PHASE 4: SKIPPED (no cell data) ---")

    # ===== PHASE 5: FINALIZE =====
    logging.info("\n--- PHASE 5: FINALIZE ---")

    # Save versions file
    logging.info("Saving versions file...")
    save_versions_file(data["outdir"])

    _timings["total"] = time.time() - t_total_start
    _timings["accounted"] = sum(v for k, v in _timings.items() if k != "total")
    _timings["unaccounted"] = _timings["total"] - _timings["accounted"]

    # Write profiling JSON
    profiling_json = data["outdir"] / "profiling_timings.json"
    with open(profiling_json, "w") as f:
        json.dump({k: round(v, 2) for k, v in _timings.items()}, f, indent=2)
    logging.info(f"Saved profiling timings to {profiling_json}")

    # Print profiling summary sorted by duration
    logging.info("\n" + "=" * 60)
    logging.info("PROFILING SUMMARY (sorted by wall-clock seconds)")
    logging.info("=" * 60)
    for k, v in sorted(
        (
            (k, v)
            for k, v in _timings.items()
            if k not in ("total", "accounted", "unaccounted")
        ),
        key=lambda x: x[1],
        reverse=True,
    ):
        pct = v / _timings["total"] * 100 if _timings["total"] > 0 else 0
        logging.info(f"  {k:40s} {v:8.1f}s  ({pct:5.1f}%)")
    logging.info(
        f"  {'--- accounted ---':40s} {_timings['accounted']:8.1f}s  ({_timings['accounted'] / _timings['total'] * 100:.1f}%)"
    )
    logging.info(
        f"  {'--- unaccounted (overhead/gaps) ---':40s} {_timings['unaccounted']:8.1f}s  ({_timings['unaccounted'] / _timings['total'] * 100:.1f}%)"
    )
    logging.info(f"  {'TOTAL':40s} {_timings['total']:8.1f}s")
    logging.info("=" * 60)

    logging.info(f"\n[TIMING] Total pipeline time: {_timings['total']:.1f}s")
    logging.info("=" * 60)
    logging.info("Combined Image QC Pipeline completed successfully!")
    logging.info(f"  Tile figures: {data['figures_dir']}")
    if has_cell_data:
        logging.info(f"  Cell figures: {data['outdir'] / 'figures'}")
        logging.info(f"  Cell metrics: {data['outdir'] / 'image_qc_metrics.json'}")
        logging.info(f"  Cell CSV: {data['outdir'] / 'image_qc_metrics.csv'}")
    logging.info(f"  Tile metrics: {data['outdir'] / 'roi_qc_metrics.json'}")
    logging.info(f"  Versions: {data['outdir'] / 'versions.yml'}")
    logging.info("=" * 60)


if __name__ == "__main__":
    main()
