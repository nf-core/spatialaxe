#!/usr/bin/env python3
"""Divide a Xenium transcripts.parquet file into spatial patches for tiled segmentation.

Standalone script — no imports from xenium_patch or any local package.
Only uses stdlib + pyarrow + numpy.

Two grid modes:
  - Uniform (default): equal-sized tiles based on --tile-width
  - Quadtree (--balanced): starts uniform, recursively subdivides dense tiles
"""

from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

XENIUM_PIXEL_SIZE_UM: float = 0.2125

TRANSCRIPT_COLS = [
    "transcript_id",
    "cell_id",
    "overlaps_nucleus",
    "feature_name",
    "x_location",
    "y_location",
    "z_location",
    "qv",
]

# Quadtree defaults
QUADTREE_MIN_TILE_WIDTH_UM: float = 200.0
QUADTREE_MAX_DEPTH: int = 4
QUADTREE_HISTOGRAM_BINS: int = 500

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Bounds:
    """Axis-aligned bounding box in either pixel or micron coordinates."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float

    @property
    def width(self) -> float:
        return self.x_max - self.x_min

    @property
    def height(self) -> float:
        return self.y_max - self.y_min


@dataclass(frozen=True)
class PatchInfo:
    """Metadata for a single patch in the grid."""

    patch_id: str
    row: int
    col: int
    global_bounds_px: Bounds
    global_bounds_um: Bounds
    core_bounds_px: Bounds
    core_bounds_um: Bounds


# ---------------------------------------------------------------------------
# Grid computation — uniform
# ---------------------------------------------------------------------------


def _compute_uniform_grid(
    image_height_px: int,
    image_width_px: int,
    grid_rows: int,
    grid_cols: int,
    overlap_px: int,
    pixel_size_um: float,
) -> list[PatchInfo]:
    """
    Compute a regular NxM grid of overlapping patches.

    Grid is computed in pixel space.  Each patch overlaps its neighbors by
    overlap_px pixels.  Core regions are computed such that every pixel
    belongs to exactly one core.

    Args:
        image_height_px: Image height in pixels.
        image_width_px: Image width in pixels.
        grid_rows: Number of rows in the patch grid.
        grid_cols: Number of columns in the patch grid.
        overlap_px: Overlap between adjacent patches in pixels.
        pixel_size_um: Microns per pixel.

    Returns:
        List of PatchInfo for every patch.
    """
    step_x = (image_width_px - overlap_px) / grid_cols
    step_y = (image_height_px - overlap_px) / grid_rows

    patches: list[PatchInfo] = []
    for row in range(grid_rows):
        for col in range(grid_cols):
            x_min_px = int(round(col * step_x))
            y_min_px = int(round(row * step_y))
            x_max_px = min(
                int(round(col * step_x + step_x + overlap_px)), image_width_px
            )
            y_max_px = min(
                int(round(row * step_y + step_y + overlap_px)), image_height_px
            )

            global_bounds_px = Bounds(x_min_px, x_max_px, y_min_px, y_max_px)

            # Core bounds: trim half-overlap from sides that have neighbors
            half_overlap = overlap_px // 2
            remainder = overlap_px % 2
            core_x_min = x_min_px + (half_overlap + remainder if col > 0 else 0)
            core_x_max = x_max_px - (half_overlap if col < grid_cols - 1 else 0)
            core_y_min = y_min_px + (half_overlap + remainder if row > 0 else 0)
            core_y_max = y_max_px - (half_overlap if row < grid_rows - 1 else 0)

            core_bounds_px = Bounds(core_x_min, core_x_max, core_y_min, core_y_max)

            global_bounds_um = Bounds(
                x_min_px * pixel_size_um,
                x_max_px * pixel_size_um,
                y_min_px * pixel_size_um,
                y_max_px * pixel_size_um,
            )
            core_bounds_um = Bounds(
                core_x_min * pixel_size_um,
                core_x_max * pixel_size_um,
                core_y_min * pixel_size_um,
                core_y_max * pixel_size_um,
            )

            patches.append(
                PatchInfo(
                    patch_id=f"patch_{row}_{col}",
                    row=row,
                    col=col,
                    global_bounds_px=global_bounds_px,
                    global_bounds_um=global_bounds_um,
                    core_bounds_px=core_bounds_px,
                    core_bounds_um=core_bounds_um,
                )
            )

    return patches


def compute_tilewidth_uniform_grid(
    image_height_px: int,
    image_width_px: int,
    tile_width_um: float,
    overlap_um: float,
    pixel_size_um: float,
    transcript_extent_um: Bounds,
) -> tuple[list[PatchInfo], int, int, int]:
    """
    Compute a uniform grid from a tile width in microns.

    Args:
        image_height_px: Image height in pixels.
        image_width_px: Image width in pixels.
        tile_width_um: Desired tile width in microns.
        overlap_um: Overlap between adjacent patches in microns.
        pixel_size_um: Size of one pixel in microns.
        transcript_extent_um: Bounding box of transcript coordinates.

    Returns:
        Tuple of (patches, grid_rows, grid_cols, overlap_px).
    """
    image_width_um = image_width_px * pixel_size_um
    image_height_um = image_height_px * pixel_size_um
    cols = max(1, math.ceil(image_width_um / tile_width_um))
    rows = max(1, math.ceil(image_height_um / tile_width_um))
    overlap_px = int(math.ceil(overlap_um / pixel_size_um))

    patches = _compute_uniform_grid(
        image_height_px, image_width_px, rows, cols, overlap_px, pixel_size_um
    )
    return patches, rows, cols, overlap_px


# ---------------------------------------------------------------------------
# Grid computation — density quadtree
# ---------------------------------------------------------------------------


def _build_prefix_sum(
    x_coords_um: np.ndarray,
    y_coords_um: np.ndarray,
    n_bins: int = QUADTREE_HISTOGRAM_BINS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build a 2D histogram and its prefix sum for fast rectangle count queries.

    Args:
        x_coords_um: Transcript X coordinates in microns.
        y_coords_um: Transcript Y coordinates in microns.
        n_bins: Number of bins along each axis.

    Returns:
        Tuple of (prefix_sum, x_edges, y_edges).
    """
    x_min, x_max = float(np.min(x_coords_um)), float(np.max(x_coords_um))
    y_min, y_max = float(np.min(y_coords_um)), float(np.max(y_coords_um))

    eps = 1e-6
    x_edges = np.linspace(x_min, x_max + eps, n_bins + 1)
    y_edges = np.linspace(y_min, y_max + eps, n_bins + 1)

    hist, _, _ = np.histogram2d(x_coords_um, y_coords_um, bins=[x_edges, y_edges])
    # hist shape is (n_bins_x, n_bins_y), transpose to (y, x) for row-major access
    hist = hist.T

    prefix_sum = np.cumsum(np.cumsum(hist, axis=0), axis=1)
    return prefix_sum, x_edges, y_edges


def _count_transcripts_in_rect(
    prefix_sum: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    x_min_um: float,
    x_max_um: float,
    y_min_um: float,
    y_max_um: float,
) -> int:
    """
    Count transcripts in a rectangle using a 2D prefix sum array.

    Args:
        prefix_sum: 2D cumulative sum array (n_bins_y x n_bins_x).
        x_edges: Histogram bin edges along X.
        y_edges: Histogram bin edges along Y.
        x_min_um: Left bound in microns.
        x_max_um: Right bound in microns.
        y_min_um: Top bound in microns.
        y_max_um: Bottom bound in microns.

    Returns:
        Approximate transcript count in the rectangle.
    """
    col_lo = max(0, int(np.searchsorted(x_edges, x_min_um, side="right")) - 1)
    col_hi = min(
        len(x_edges) - 1, int(np.searchsorted(x_edges, x_max_um, side="right")) - 1
    )
    row_lo = max(0, int(np.searchsorted(y_edges, y_min_um, side="right")) - 1)
    row_hi = min(
        len(y_edges) - 1, int(np.searchsorted(y_edges, y_max_um, side="right")) - 1
    )

    col_hi = min(col_hi, prefix_sum.shape[1] - 1)
    row_hi = min(row_hi, prefix_sum.shape[0] - 1)

    if col_lo > col_hi or row_lo > row_hi:
        return 0

    total = int(
        prefix_sum[row_hi, col_hi]
        - (prefix_sum[row_lo - 1, col_hi] if row_lo > 0 else 0)
        - (prefix_sum[row_hi, col_lo - 1] if col_lo > 0 else 0)
        + (prefix_sum[row_lo - 1, col_lo - 1] if row_lo > 0 and col_lo > 0 else 0)
    )
    return max(0, total)


def _subdivide_regions(
    regions: list[tuple[float, float, float, float]],
    prefix_sum: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    max_transcripts: int,
    min_tile_width_um: float,
    max_depth: int,
) -> list[tuple[float, float, float, float]]:
    """
    Recursively subdivide regions exceeding the transcript threshold.

    Uses a stack instead of recursion for large grids.

    Args:
        regions: List of (x_min, x_max, y_min, y_max) tuples in microns.
        prefix_sum: 2D prefix sum for fast counting.
        x_edges: Histogram X bin edges.
        y_edges: Histogram Y bin edges.
        max_transcripts: Maximum transcripts allowed per region.
        min_tile_width_um: Minimum tile dimension before stopping.
        max_depth: Maximum recursion depth.

    Returns:
        List of final (x_min, x_max, y_min, y_max) regions.
    """
    result: list[tuple[float, float, float, float]] = []
    stack: list[tuple[tuple[float, float, float, float], int]] = [
        (r, 0) for r in regions
    ]

    while stack:
        region, depth = stack.pop()
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min
        height = y_max - y_min

        count = _count_transcripts_in_rect(
            prefix_sum, x_edges, y_edges, x_min, x_max, y_min, y_max
        )

        if count <= max_transcripts or depth >= max_depth:
            result.append(region)
            continue

        if min(width, height) / 2 < min_tile_width_um:
            result.append(region)
            continue

        # Split into 4 quadrants
        mid_x = (x_min + x_max) / 2
        mid_y = (y_min + y_max) / 2
        children = [
            (x_min, mid_x, y_min, mid_y),
            (mid_x, x_max, y_min, mid_y),
            (x_min, mid_x, mid_y, y_max),
            (mid_x, x_max, mid_y, y_max),
        ]
        for child in children:
            stack.append((child, depth + 1))

    return result


def _regions_to_patches(
    regions: list[tuple[float, float, float, float]],
    overlap_um: float,
    overlap_px: int,
    pixel_size_um: float,
    image_width_px: int,
    image_height_px: int,
) -> list[PatchInfo]:
    """
    Convert quadtree regions to PatchInfo objects with overlap.

    Args:
        regions: Sorted list of (x_min, x_max, y_min, y_max) in microns.
        overlap_um: Overlap in microns.
        overlap_px: Overlap in pixels.
        pixel_size_um: Microns per pixel.
        image_width_px: Image width in pixels.
        image_height_px: Image height in pixels.

    Returns:
        List of PatchInfo objects.
    """
    patches: list[PatchInfo] = []
    for i, (x_min_um, x_max_um, y_min_um, y_max_um) in enumerate(regions):
        # Core bounds in pixels
        core_x_min_px = max(
            0, min(int(round(x_min_um / pixel_size_um)), image_width_px)
        )
        core_x_max_px = max(
            0, min(int(round(x_max_um / pixel_size_um)), image_width_px)
        )
        core_y_min_px = max(
            0, min(int(round(y_min_um / pixel_size_um)), image_height_px)
        )
        core_y_max_px = max(
            0, min(int(round(y_max_um / pixel_size_um)), image_height_px)
        )

        core_bounds_px = Bounds(
            core_x_min_px, core_x_max_px, core_y_min_px, core_y_max_px
        )

        # Global bounds: core extended by overlap, clamped to image
        global_x_min_px = max(0, core_x_min_px - overlap_px)
        global_x_max_px = min(image_width_px, core_x_max_px + overlap_px)
        global_y_min_px = max(0, core_y_min_px - overlap_px)
        global_y_max_px = min(image_height_px, core_y_max_px + overlap_px)

        global_bounds_px = Bounds(
            global_x_min_px, global_x_max_px, global_y_min_px, global_y_max_px
        )

        core_bounds_um = Bounds(
            core_x_min_px * pixel_size_um,
            core_x_max_px * pixel_size_um,
            core_y_min_px * pixel_size_um,
            core_y_max_px * pixel_size_um,
        )
        global_bounds_um = Bounds(
            global_x_min_px * pixel_size_um,
            global_x_max_px * pixel_size_um,
            global_y_min_px * pixel_size_um,
            global_y_max_px * pixel_size_um,
        )

        patches.append(
            PatchInfo(
                patch_id=f"patch_{i}",
                row=i,
                col=0,
                global_bounds_px=global_bounds_px,
                global_bounds_um=global_bounds_um,
                core_bounds_px=core_bounds_px,
                core_bounds_um=core_bounds_um,
            )
        )

    return patches


def compute_density_quadtree_grid(
    image_height_px: int,
    image_width_px: int,
    tile_width_um: float,
    overlap_um: float,
    pixel_size_um: float,
    x_coords_um: np.ndarray,
    y_coords_um: np.ndarray,
    max_transcripts_per_patch: int | None = None,
    min_tile_width_um: float = QUADTREE_MIN_TILE_WIDTH_UM,
    max_depth: int = QUADTREE_MAX_DEPTH,
) -> tuple[list[PatchInfo], int, int, int]:
    """
    Compute an adaptive quadtree grid that subdivides dense regions.

    Starts with a uniform grid derived from tile_width_um, then recursively
    subdivides patches exceeding max_transcripts_per_patch.

    Args:
        image_height_px: Image height in pixels.
        image_width_px: Image width in pixels.
        tile_width_um: Base tile width in microns.
        overlap_um: Overlap between adjacent patches in microns.
        pixel_size_um: Microns per pixel.
        x_coords_um: Transcript X coordinates in microns.
        y_coords_um: Transcript Y coordinates in microns.
        max_transcripts_per_patch: Target max transcripts per patch.
            If None, auto-computed as 2x the average per initial patch.
        min_tile_width_um: Minimum tile dimension before stopping.
        max_depth: Maximum recursion depth.

    Returns:
        Tuple of (patches, initial_rows, initial_cols, overlap_px).
    """
    image_width_um = image_width_px * pixel_size_um
    image_height_um = image_height_px * pixel_size_um
    overlap_px = int(math.ceil(overlap_um / pixel_size_um))

    initial_cols = max(1, math.ceil(image_width_um / tile_width_um))
    initial_rows = max(1, math.ceil(image_height_um / tile_width_um))

    # Build prefix sum for fast counting
    prefix_sum, x_edges, y_edges = _build_prefix_sum(x_coords_um, y_coords_um)

    # Define initial regions in microns
    cell_width_um = image_width_um / initial_cols
    cell_height_um = image_height_um / initial_rows

    initial_regions: list[tuple[float, float, float, float]] = []
    for row in range(initial_rows):
        for col in range(initial_cols):
            x_min = col * cell_width_um
            x_max = min((col + 1) * cell_width_um, image_width_um)
            y_min = row * cell_height_um
            y_max = min((row + 1) * cell_height_um, image_height_um)
            initial_regions.append((x_min, x_max, y_min, y_max))

    # Auto-compute threshold
    n_initial = len(initial_regions)
    total_transcripts = len(x_coords_um)
    if max_transcripts_per_patch is None:
        max_transcripts_per_patch = max(1, int(total_transcripts / n_initial * 2))

    # Recursive subdivision
    final_regions = _subdivide_regions(
        initial_regions,
        prefix_sum,
        x_edges,
        y_edges,
        max_transcripts_per_patch,
        min_tile_width_um,
        max_depth,
    )

    # Sort by (y_min, x_min) for deterministic ordering
    final_regions.sort(key=lambda r: (r[2], r[0]))

    # Convert to PatchInfo
    patches = _regions_to_patches(
        final_regions,
        overlap_um,
        overlap_px,
        pixel_size_um,
        image_width_px,
        image_height_px,
    )

    return patches, initial_rows, initial_cols, overlap_px


# ---------------------------------------------------------------------------
# Sparse tile merging
# ---------------------------------------------------------------------------


def _count_transcripts_per_tile(
    patches: list[PatchInfo],
    x_coords_um: np.ndarray,
    y_coords_um: np.ndarray,
) -> dict[str, int]:
    """
    Count transcripts falling within each patch's core bounds.

    Uses core bounds (not global) to avoid double-counting transcripts
    in overlap regions.

    Args:
        patches: List of PatchInfo objects.
        x_coords_um: Transcript X coordinates in microns.
        y_coords_um: Transcript Y coordinates in microns.

    Returns:
        Dict mapping patch_id to transcript count.
    """
    counts: dict[str, int] = {}
    for p in patches:
        cb = p.core_bounds_um
        mask = (
            (x_coords_um >= cb.x_min)
            & (x_coords_um < cb.x_max)
            & (y_coords_um >= cb.y_min)
            & (y_coords_um < cb.y_max)
        )
        counts[p.patch_id] = int(np.sum(mask))
    return counts


def _find_adjacent_patches(
    patches: list[PatchInfo],
) -> dict[str, list[str]]:
    """
    Build an adjacency map: patches sharing a core bounds edge are neighbors.

    Two patches are adjacent if their core bounds share an edge (touch or
    overlap along one axis while overlapping along the other axis).

    Args:
        patches: List of PatchInfo objects.

    Returns:
        Dict mapping patch_id to list of adjacent patch_ids.
    """
    adjacency: dict[str, list[str]] = {p.patch_id: [] for p in patches}
    eps = 1.0  # tolerance in microns for edge sharing

    for i, a in enumerate(patches):
        for j in range(i + 1, len(patches)):
            b = patches[j]
            ac = a.core_bounds_um
            bc = b.core_bounds_um

            # Check X-axis overlap (cores overlap in X)
            x_overlap = ac.x_min < bc.x_max and bc.x_min < ac.x_max
            # Check Y-axis overlap (cores overlap in Y)
            y_overlap = ac.y_min < bc.y_max and bc.y_min < ac.y_max

            # Adjacent along X: share a vertical edge, overlap in Y
            x_touching = (
                abs(ac.x_max - bc.x_min) < eps or abs(bc.x_max - ac.x_min) < eps
            )
            # Adjacent along Y: share a horizontal edge, overlap in X
            y_touching = (
                abs(ac.y_max - bc.y_min) < eps or abs(bc.y_max - ac.y_min) < eps
            )

            if (x_touching and y_overlap) or (y_touching and x_overlap):
                adjacency[a.patch_id].append(b.patch_id)
                adjacency[b.patch_id].append(a.patch_id)

    return adjacency


def _recalculate_core_bounds(
    patches: list[PatchInfo],
    overlap_px: int,
    pixel_size_um: float,
    image_width_px: int,
    image_height_px: int,
) -> list[PatchInfo]:
    """
    Recalculate core bounds for all patches after merging.

    Core bounds are derived from the regions: the core is the
    non-overlapping portion of each tile. After merging, we extract
    core regions from global bounds by trimming the overlap, then
    rebuild PatchInfo objects.

    For merged grids where tiles may be irregular, core bounds equal
    the global bounds shrunk by half the overlap on each side that has
    a neighbor, clamped to the image extent.

    Args:
        patches: Current list of PatchInfo (with updated global bounds).
        overlap_px: Overlap in pixels.
        pixel_size_um: Microns per pixel.
        image_width_px: Image width in pixels.
        image_height_px: Image height in pixels.

    Returns:
        New list of PatchInfo with recalculated core and global bounds.
    """
    if not patches:
        return []

    # Extract core regions in microns from global bounds minus overlap
    half_overlap_um = (overlap_px * pixel_size_um) / 2.0
    image_width_um = image_width_px * pixel_size_um
    image_height_um = image_height_px * pixel_size_um

    # Collect all core regions (global shrunk by half overlap)
    core_regions_um: list[tuple[float, float, float, float]] = []
    for p in patches:
        gb = p.global_bounds_um
        # Shrink by half overlap on each side, but not past image edge
        cx_min = gb.x_min + (half_overlap_um if gb.x_min > 0 else 0)
        cx_max = gb.x_max - (half_overlap_um if gb.x_max < image_width_um else 0)
        cy_min = gb.y_min + (half_overlap_um if gb.y_min > 0 else 0)
        cy_max = gb.y_max - (half_overlap_um if gb.y_max < image_height_um else 0)
        core_regions_um.append((cx_min, cx_max, cy_min, cy_max))

    # Rebuild patches using core regions -> global bounds (core + overlap)
    result: list[PatchInfo] = []
    for i, p in enumerate(patches):
        cx_min, cx_max, cy_min, cy_max = core_regions_um[i]

        # Core bounds in pixels
        core_x_min_px = max(0, min(int(round(cx_min / pixel_size_um)), image_width_px))
        core_x_max_px = max(0, min(int(round(cx_max / pixel_size_um)), image_width_px))
        core_y_min_px = max(0, min(int(round(cy_min / pixel_size_um)), image_height_px))
        core_y_max_px = max(0, min(int(round(cy_max / pixel_size_um)), image_height_px))

        core_bounds_px = Bounds(
            core_x_min_px, core_x_max_px, core_y_min_px, core_y_max_px
        )

        # Global bounds: core extended by overlap, clamped to image
        global_x_min_px = max(0, core_x_min_px - overlap_px)
        global_x_max_px = min(image_width_px, core_x_max_px + overlap_px)
        global_y_min_px = max(0, core_y_min_px - overlap_px)
        global_y_max_px = min(image_height_px, core_y_max_px + overlap_px)

        global_bounds_px = Bounds(
            global_x_min_px, global_x_max_px, global_y_min_px, global_y_max_px
        )

        core_bounds_um = Bounds(
            core_x_min_px * pixel_size_um,
            core_x_max_px * pixel_size_um,
            core_y_min_px * pixel_size_um,
            core_y_max_px * pixel_size_um,
        )
        global_bounds_um = Bounds(
            global_x_min_px * pixel_size_um,
            global_x_max_px * pixel_size_um,
            global_y_min_px * pixel_size_um,
            global_y_max_px * pixel_size_um,
        )

        result.append(
            PatchInfo(
                patch_id=p.patch_id,
                row=p.row,
                col=p.col,
                global_bounds_px=global_bounds_px,
                global_bounds_um=global_bounds_um,
                core_bounds_px=core_bounds_px,
                core_bounds_um=core_bounds_um,
            )
        )

    return result


def merge_sparse_tiles(
    patches: list[PatchInfo],
    x_coords_um: np.ndarray,
    y_coords_um: np.ndarray,
    overlap_px: int,
    pixel_size_um: float,
    image_width_px: int,
    image_height_px: int,
    min_transcripts: int = 1000,
) -> tuple[list[PatchInfo], int]:
    """
    Merge tiles below min_transcripts into their least populated adjacent neighbor.

    Iteratively finds the sparsest tile below the threshold and merges it
    into its smallest neighbor for balanced tile sizes. Repeats until no
    tiles remain below the threshold (or a tile has no neighbors to merge into).

    Args:
        patches: List of PatchInfo objects from grid computation.
        x_coords_um: Transcript X coordinates in microns.
        y_coords_um: Transcript Y coordinates in microns.
        overlap_px: Overlap in pixels.
        pixel_size_um: Microns per pixel.
        image_width_px: Image width in pixels.
        image_height_px: Image height in pixels.
        min_transcripts: Minimum transcript count per tile.

    Returns:
        Tuple of (merged patches, number of merges performed).
    """
    if len(patches) <= 1:
        return patches, 0

    # Work with mutable list
    active = list(patches)
    merge_count = 0

    while True:
        counts = _count_transcripts_per_tile(active, x_coords_um, y_coords_um)
        adjacency = _find_adjacent_patches(active)

        # Find sparsest tile below threshold
        sparse_candidates = [
            (pid, cnt) for pid, cnt in counts.items() if cnt < min_transcripts
        ]
        if not sparse_candidates:
            break

        # Sort by count ascending to merge sparsest first
        sparse_candidates.sort(key=lambda t: t[1])
        sparse_id, sparse_count = sparse_candidates[0]

        # Find neighbors and pick the least populated one for balanced merging
        neighbors = adjacency.get(sparse_id, [])
        if not neighbors:
            # No neighbors for this tile — skip it and try next sparsest
            sparse_candidates = [(pid, cnt) for pid, cnt in sparse_candidates[1:]]
            found = False
            for pid, cnt in sparse_candidates:
                nbrs = adjacency.get(pid, [])
                if nbrs:
                    sparse_id, sparse_count = pid, cnt
                    neighbors = nbrs
                    found = True
                    break
            if not found:
                break

        best_neighbor_id = min(neighbors, key=lambda nid: counts.get(nid, 0))

        # Find the actual PatchInfo objects
        sparse_patch = next(p for p in active if p.patch_id == sparse_id)
        neighbor_patch = next(p for p in active if p.patch_id == best_neighbor_id)

        # Expand neighbor's global bounds to cover both tiles
        sg = sparse_patch.global_bounds_um
        ng = neighbor_patch.global_bounds_um
        merged_global_um = Bounds(
            x_min=min(sg.x_min, ng.x_min),
            x_max=max(sg.x_max, ng.x_max),
            y_min=min(sg.y_min, ng.y_min),
            y_max=max(sg.y_max, ng.y_max),
        )

        # Also merge core bounds (union)
        sc = sparse_patch.core_bounds_um
        nc = neighbor_patch.core_bounds_um
        merged_core_um = Bounds(
            x_min=min(sc.x_min, nc.x_min),
            x_max=max(sc.x_max, nc.x_max),
            y_min=min(sc.y_min, nc.y_min),
            y_max=max(sc.y_max, nc.y_max),
        )

        # Convert merged bounds to pixels
        merged_global_px = Bounds(
            x_min=max(0, int(round(merged_global_um.x_min / pixel_size_um))),
            x_max=min(
                image_width_px, int(round(merged_global_um.x_max / pixel_size_um))
            ),
            y_min=max(0, int(round(merged_global_um.y_min / pixel_size_um))),
            y_max=min(
                image_height_px, int(round(merged_global_um.y_max / pixel_size_um))
            ),
        )
        merged_core_px = Bounds(
            x_min=max(0, int(round(merged_core_um.x_min / pixel_size_um))),
            x_max=min(image_width_px, int(round(merged_core_um.x_max / pixel_size_um))),
            y_min=max(0, int(round(merged_core_um.y_min / pixel_size_um))),
            y_max=min(
                image_height_px, int(round(merged_core_um.y_max / pixel_size_um))
            ),
        )

        # Create merged patch (keeps absorbing tile's ID and position)
        merged_patch = PatchInfo(
            patch_id=neighbor_patch.patch_id,
            row=neighbor_patch.row,
            col=neighbor_patch.col,
            global_bounds_px=merged_global_px,
            global_bounds_um=merged_global_um,
            core_bounds_px=merged_core_px,
            core_bounds_um=merged_core_um,
        )

        # Replace neighbor with merged patch and remove sparse tile
        active = [
            merged_patch if p.patch_id == best_neighbor_id else p
            for p in active
            if p.patch_id != sparse_id
        ]
        merge_count += 1

        print(
            f"  Merged {sparse_id} ({sparse_count:,} transcripts) "
            f"into {best_neighbor_id} ({counts[best_neighbor_id]:,} transcripts)"
        )

    if merge_count > 0:
        # Recalculate core bounds for consistency
        active = _recalculate_core_bounds(
            active, overlap_px, pixel_size_um, image_width_px, image_height_px
        )

    return active, merge_count


# ---------------------------------------------------------------------------
# Transcript division
# ---------------------------------------------------------------------------


def _filter_and_write_patch_transcripts(
    full_table: pa.Table,
    output_path: Path,
    bounds_um: Bounds,
    origin_x: float,
    origin_y: float,
) -> int:
    """
    Filter transcripts to a spatial region and write to parquet.

    Transcripts are filtered to global_bounds (including overlap), then
    coordinates are offset by subtracting the global_bounds origin.

    Args:
        full_table: Full transcript table as a pyarrow Table.
        output_path: Path for the filtered output parquet.
        bounds_um: Spatial bounding box for filtering (microns).
        origin_x: X offset to subtract for local coordinates.
        origin_y: Y offset to subtract for local coordinates.

    Returns:
        Number of transcripts written.
    """
    x_col = full_table.column("x_location")
    y_col = full_table.column("y_location")

    mask = pc.and_(
        pc.and_(
            pc.greater_equal(x_col, pa.scalar(bounds_um.x_min, type=x_col.type)),
            pc.less(x_col, pa.scalar(bounds_um.x_max, type=x_col.type)),
        ),
        pc.and_(
            pc.greater_equal(y_col, pa.scalar(bounds_um.y_min, type=y_col.type)),
            pc.less(y_col, pa.scalar(bounds_um.y_max, type=y_col.type)),
        ),
    )
    filtered = full_table.filter(mask)

    if origin_x != 0.0 or origin_y != 0.0:
        fx = filtered.column("x_location")
        fy = filtered.column("y_location")
        x_local = pc.subtract(fx, pa.scalar(origin_x, type=fx.type))
        y_local = pc.subtract(fy, pa.scalar(origin_y, type=fy.type))
        idx_x = filtered.schema.get_field_index("x_location")
        idx_y = filtered.schema.get_field_index("y_location")
        filtered = filtered.set_column(idx_x, "x_location", x_local)
        filtered = filtered.set_column(idx_y, "y_location", y_local)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(filtered, str(output_path))
    return len(filtered)


def _process_patch(
    patch: PatchInfo,
    output_dir: Path,
    full_table: pa.Table,
) -> int:
    """
    Write transcript subset for a single patch.

    Args:
        patch: Patch metadata.
        output_dir: Root output directory.
        full_table: Full transcript table.

    Returns:
        Number of transcripts written.
    """
    patch_dir = output_dir / patch.patch_id
    bounds_um = patch.global_bounds_um
    return _filter_and_write_patch_transcripts(
        full_table,
        patch_dir / "transcripts.parquet",
        bounds_um,
        origin_x=bounds_um.x_min,
        origin_y=bounds_um.y_min,
    )


# ---------------------------------------------------------------------------
# JSON serialization
# ---------------------------------------------------------------------------


def _bounds_to_dict(b: Bounds) -> dict[str, float]:
    """Serialize a Bounds to a JSON-compatible dict."""
    return {"x_min": b.x_min, "x_max": b.x_max, "y_min": b.y_min, "y_max": b.y_max}


def save_grid_metadata(
    patches: list[PatchInfo],
    image_height_px: int,
    image_width_px: int,
    pixel_size_um: float,
    transcript_extent_um: Bounds,
    grid_rows: int,
    grid_cols: int,
    overlap_um: float,
    overlap_px: int,
    grid_type: str,
    output_path: Path,
) -> None:
    """
    Serialize grid metadata to JSON.

    Args:
        patches: List of PatchInfo objects.
        image_height_px: Image height in pixels.
        image_width_px: Image width in pixels.
        pixel_size_um: Microns per pixel.
        transcript_extent_um: Bounding box of transcript coordinates.
        grid_rows: Number of rows in the initial grid.
        grid_cols: Number of columns in the initial grid.
        overlap_um: Overlap in microns.
        overlap_px: Overlap in pixels.
        grid_type: Grid type string ("uniform" or "density_quadtree").
        output_path: Path to write JSON file.
    """
    data = {
        "version": "1.0",
        "bundle_path": "",
        "image_height_px": image_height_px,
        "image_width_px": image_width_px,
        "pixel_size_um": pixel_size_um,
        "transcript_extent_um": _bounds_to_dict(transcript_extent_um),
        "grid_rows": grid_rows,
        "grid_cols": grid_cols,
        "overlap_um": overlap_um,
        "overlap_px": overlap_px,
        "grid_type": grid_type,
        "patches": [
            {
                "patch_id": p.patch_id,
                "row": p.row,
                "col": p.col,
                "global_bounds_px": _bounds_to_dict(p.global_bounds_px),
                "global_bounds_um": _bounds_to_dict(p.global_bounds_um),
                "core_bounds_px": _bounds_to_dict(p.core_bounds_px),
                "core_bounds_um": _bounds_to_dict(p.core_bounds_um),
            }
            for p in patches
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


# ---------------------------------------------------------------------------
# Coordinate shift helper
# ---------------------------------------------------------------------------


def _shift_patches_to_real_coords(
    patches: list[PatchInfo],
    ox: float,
    oy: float,
) -> list[PatchInfo]:
    """
    Shift patch micron bounds by (ox, oy) to align with real transcript coords.

    Pixel bounds remain zero-origin (there is no real image to index into).

    Args:
        patches: Patches in zero-origin micron space.
        ox: X offset (transcript extent x_min).
        oy: Y offset (transcript extent y_min).

    Returns:
        New list of PatchInfo with shifted micron bounds.
    """
    shifted: list[PatchInfo] = []
    for p in patches:
        gu = p.global_bounds_um
        cu = p.core_bounds_um
        shifted.append(
            PatchInfo(
                patch_id=p.patch_id,
                row=p.row,
                col=p.col,
                global_bounds_px=p.global_bounds_px,
                global_bounds_um=Bounds(
                    gu.x_min + ox, gu.x_max + ox, gu.y_min + oy, gu.y_max + oy
                ),
                core_bounds_px=p.core_bounds_px,
                core_bounds_um=Bounds(
                    cu.x_min + ox, cu.x_max + ox, cu.y_min + oy, cu.y_max + oy
                ),
            )
        )
    return shifted


# ---------------------------------------------------------------------------
# Main divide logic
# ---------------------------------------------------------------------------


def divide_transcripts(
    transcripts_path: Path,
    output_dir: Path,
    image_width_px: int,
    image_height_px: int,
    tile_width_um: float,
    overlap_um: float,
    balanced: bool,
    pixel_size_um: float = XENIUM_PIXEL_SIZE_UM,
    max_workers: int | None = None,
    min_transcripts: int = 1000,
) -> None:
    """
    Divide transcripts into overlapping spatial patches.

    Reads the transcript table once, computes a grid, merges sparse tiles
    into neighbors, and writes per-patch parquet files with coordinates
    offset to patch-local space.

    Args:
        transcripts_path: Path to transcripts.parquet.
        output_dir: Output directory for patches.
        image_width_px: Image width in pixels.
        image_height_px: Image height in pixels.
        tile_width_um: Tile width in microns.
        overlap_um: Overlap between adjacent patches in microns.
        balanced: If True, use density quadtree mode.
        pixel_size_um: Microns per pixel.
        max_workers: Maximum threads for parallel patch writes.
        min_transcripts: Minimum transcripts per tile; sparse tiles merged
            into neighbors. Set to 0 to disable merging.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read full transcript table
    full_table = pq.read_table(str(transcripts_path))
    n_total = len(full_table)
    print(f"Read {n_total:,} transcripts from {transcripts_path}")

    # Compute transcript extent
    x_col = full_table.column("x_location")
    y_col = full_table.column("y_location")
    extent_um = Bounds(
        x_min=pc.min(x_col).as_py(),
        x_max=pc.max(x_col).as_py(),
        y_min=pc.min(y_col).as_py(),
        y_max=pc.max(y_col).as_py(),
    )
    print(
        f"Transcript extent: "
        f"x=[{extent_um.x_min:.1f}, {extent_um.x_max:.1f}] "
        f"y=[{extent_um.y_min:.1f}, {extent_um.y_max:.1f}] um"
    )

    # Build grid in zero-origin space when transcripts have a positive offset.
    # The grid functions work in pixel space starting at (0, 0).  We shift
    # micron bounds back to real coordinates afterward.
    ox = extent_um.x_min
    oy = extent_um.y_min

    if balanced:
        # Shift coordinates to zero-origin for density computation
        x_coords = x_col.to_numpy() - ox
        y_coords = y_col.to_numpy() - oy

        patches, grid_rows, grid_cols, overlap_px = compute_density_quadtree_grid(
            image_height_px=image_height_px,
            image_width_px=image_width_px,
            tile_width_um=tile_width_um,
            overlap_um=overlap_um,
            pixel_size_um=pixel_size_um,
            x_coords_um=x_coords,
            y_coords_um=y_coords,
        )
        grid_type = "density_quadtree"
    else:
        patches, grid_rows, grid_cols, overlap_px = compute_tilewidth_uniform_grid(
            image_height_px=image_height_px,
            image_width_px=image_width_px,
            tile_width_um=tile_width_um,
            overlap_um=overlap_um,
            pixel_size_um=pixel_size_um,
            transcript_extent_um=extent_um,
        )
        grid_type = "uniform"

    # Merge sparse tiles into neighbors
    n_before_merge = len(patches)
    if min_transcripts > 0 and len(patches) > 1:
        # Coordinates for counting: use zero-origin if not already
        if balanced:
            merge_x = x_coords
            merge_y = y_coords
        else:
            merge_x = x_col.to_numpy() - ox
            merge_y = y_col.to_numpy() - oy

        patches, n_merged = merge_sparse_tiles(
            patches=patches,
            x_coords_um=merge_x,
            y_coords_um=merge_y,
            overlap_px=overlap_px,
            pixel_size_um=pixel_size_um,
            image_width_px=image_width_px,
            image_height_px=image_height_px,
            min_transcripts=min_transcripts,
        )
        if n_merged > 0:
            grid_type = f"{grid_type}+merged"
            print(
                f"Merged {n_merged} sparse tiles: "
                f"{n_before_merge} -> {len(patches)} patches"
            )

    # Shift micron bounds to real transcript coordinates
    if ox != 0.0 or oy != 0.0:
        patches = _shift_patches_to_real_coords(patches, ox, oy)

    print(
        f"Grid: {grid_type}, {grid_rows}x{grid_cols} initial, "
        f"{len(patches)} patches, overlap={overlap_um} um"
    )

    # Write patches in parallel
    n_patches = len(patches)
    workers = (
        max_workers if max_workers is not None else min(n_patches, os.cpu_count() or 1)
    )

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(_process_patch, patch, output_dir, full_table)
            for patch in patches
        ]
        for i, future in enumerate(futures):
            count = future.result()
            print(f"  {patches[i].patch_id}: {count:,} transcripts")

    # Save grid metadata
    save_grid_metadata(
        patches=patches,
        image_height_px=image_height_px,
        image_width_px=image_width_px,
        pixel_size_um=pixel_size_um,
        transcript_extent_um=extent_um,
        grid_rows=grid_rows,
        grid_cols=grid_cols,
        overlap_um=overlap_um,
        overlap_px=overlap_px,
        grid_type=grid_type,
        output_path=output_dir / "patch_grid.json",
    )
    print(f"Grid metadata saved to {output_dir / 'patch_grid.json'}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Args:
        argv: Argument list (defaults to sys.argv[1:]).

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Divide Xenium transcripts.parquet into spatial patches.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--transcripts",
        type=Path,
        required=True,
        help="Path to transcripts.parquet",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for patches",
    )
    parser.add_argument(
        "--tile-width",
        type=float,
        default=2000.0,
        help="Tile width in microns",
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=50.0,
        help="Overlap between patches in microns",
    )
    parser.add_argument(
        "--balanced",
        action="store_true",
        help="Enable density quadtree mode (subdivides dense tiles)",
    )
    parser.add_argument(
        "--image-width",
        type=int,
        required=True,
        help="Image width in pixels",
    )
    parser.add_argument(
        "--image-height",
        type=int,
        required=True,
        help="Image height in pixels",
    )
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=XENIUM_PIXEL_SIZE_UM,
        help="Pixel size in microns",
    )
    parser.add_argument(
        "--min-transcripts",
        type=int,
        default=1000,
        help="Minimum transcripts per tile; sparse tiles are merged into neighbors",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Maximum threads for parallel writes",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Entry point."""
    args = parse_args(argv)

    divide_transcripts(
        transcripts_path=args.transcripts,
        output_dir=args.output,
        image_width_px=args.image_width,
        image_height_px=args.image_height,
        tile_width_um=args.tile_width,
        overlap_um=args.overlap,
        balanced=args.balanced,
        pixel_size_um=args.pixel_size,
        max_workers=args.max_workers,
        min_transcripts=args.min_transcripts,
    )


if __name__ == "__main__":
    main()
