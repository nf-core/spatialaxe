"""Tests for divide_transcripts.py — grid computation + transcript division."""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

# ---------------------------------------------------------------------------
# Import the standalone script from bin/
# ---------------------------------------------------------------------------

_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
_SCRIPT = _BIN_DIR / "divide_transcripts.py"
_spec = importlib.util.spec_from_file_location("divide_transcripts", _SCRIPT)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["divide_transcripts"] = _mod
_spec.loader.exec_module(_mod)

from divide_transcripts import (  # noqa: E402
    Bounds,
    PatchInfo,
    _compute_uniform_grid,
    _count_transcripts_per_tile,
    _find_adjacent_patches,
    compute_density_quadtree_grid,
    compute_tilewidth_uniform_grid,
    divide_transcripts,
    merge_sparse_tiles,
    save_grid_metadata,
)

PIXEL_SIZE = 0.2125


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def synthetic_transcripts(tmp_path: Path) -> Path:
    """Write a synthetic transcripts.parquet with 1000 rows, uniform spatial distribution."""
    rng = np.random.default_rng(42)
    n = 1000
    table = pa.table(
        {
            "transcript_id": pa.array([f"tx_{i}" for i in range(n)], type=pa.string()),
            "cell_id": pa.array(["UNASSIGNED"] * n, type=pa.string()),
            "overlaps_nucleus": pa.array([0] * n, type=pa.int32()),
            "feature_name": pa.array(
                [f"gene_{i % 50}" for i in range(n)], type=pa.string()
            ),
            "x_location": pa.array(rng.uniform(0.0, 1275.0, n), type=pa.float32()),
            "y_location": pa.array(rng.uniform(0.0, 1275.0, n), type=pa.float32()),
            "z_location": pa.array(rng.uniform(0.0, 10.0, n), type=pa.float32()),
            "qv": pa.array(rng.uniform(20.0, 40.0, n), type=pa.float32()),
        }
    )
    path = tmp_path / "transcripts.parquet"
    pq.write_table(table, str(path))
    return path


@pytest.fixture
def dense_corner_transcripts(tmp_path: Path) -> Path:
    """90% of transcripts in the top-left corner, 10% uniform."""
    rng = np.random.default_rng(99)
    n_dense = 900
    n_sparse = 100
    n = n_dense + n_sparse

    x_dense = rng.uniform(0.0, 50.0, n_dense)
    y_dense = rng.uniform(0.0, 50.0, n_dense)
    x_sparse = rng.uniform(0.0, 1275.0, n_sparse)
    y_sparse = rng.uniform(0.0, 1275.0, n_sparse)

    table = pa.table(
        {
            "transcript_id": pa.array([f"tx_{i}" for i in range(n)], type=pa.string()),
            "cell_id": pa.array(["UNASSIGNED"] * n, type=pa.string()),
            "overlaps_nucleus": pa.array([0] * n, type=pa.int32()),
            "feature_name": pa.array(
                [f"gene_{i % 50}" for i in range(n)], type=pa.string()
            ),
            "x_location": pa.array(
                np.concatenate([x_dense, x_sparse]).astype(np.float32),
                type=pa.float32(),
            ),
            "y_location": pa.array(
                np.concatenate([y_dense, y_sparse]).astype(np.float32),
                type=pa.float32(),
            ),
            "z_location": pa.array(
                rng.uniform(0.0, 10.0, n).astype(np.float32), type=pa.float32()
            ),
            "qv": pa.array(
                rng.uniform(20.0, 40.0, n).astype(np.float32), type=pa.float32()
            ),
        }
    )
    path = tmp_path / "transcripts_dense.parquet"
    pq.write_table(table, str(path))
    return path


# ---------------------------------------------------------------------------
# Uniform grid tests
# ---------------------------------------------------------------------------


class TestUniformGridBasic:
    def test_uniform_grid_basic(self):
        """3x3 grid from 2000um tile width on a 6000x6000um image."""
        # Use exact pixel count: 6000um / 0.2125um/px = 28235.29... -> round to 28235
        # image_um = 28235 * 0.2125 = 5999.9375; ceil(5999.9375 / 2000) = 3
        image_px = int(6000.0 / PIXEL_SIZE)
        extent = Bounds(0.0, 6000.0, 0.0, 6000.0)

        patches, rows, cols, overlap_px = compute_tilewidth_uniform_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=2000.0,
            overlap_um=50.0,
            pixel_size_um=PIXEL_SIZE,
            transcript_extent_um=extent,
        )

        assert rows == 3
        assert cols == 3
        assert len(patches) == 9

    def test_uniform_grid_single_tile(self):
        """Tile width larger than image produces a 1x1 grid."""
        image_px = 1000
        extent = Bounds(0.0, image_px * PIXEL_SIZE, 0.0, image_px * PIXEL_SIZE)

        patches, rows, cols, _ = compute_tilewidth_uniform_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=50000.0,
            overlap_um=50.0,
            pixel_size_um=PIXEL_SIZE,
            transcript_extent_um=extent,
        )

        assert rows == 1
        assert cols == 1
        assert len(patches) == 1

    def test_uniform_grid_overlap(self):
        """Global bounds extend beyond core by the overlap amount."""
        image_px = 1000
        overlap_um = 50.0

        patches = _compute_uniform_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            grid_rows=2,
            grid_cols=2,
            overlap_px=int(math.ceil(overlap_um / PIXEL_SIZE)),
            pixel_size_um=PIXEL_SIZE,
        )

        # Interior patch boundary: global should extend beyond core
        for p in patches:
            assert p.global_bounds_px.x_min <= p.core_bounds_px.x_min
            assert p.global_bounds_px.x_max >= p.core_bounds_px.x_max
            assert p.global_bounds_px.y_min <= p.core_bounds_px.y_min
            assert p.global_bounds_px.y_max >= p.core_bounds_px.y_max


# ---------------------------------------------------------------------------
# Quadtree grid tests
# ---------------------------------------------------------------------------


class TestQuadtreeGrid:
    def test_quadtree_uniform_density(self):
        """When density is uniform with high threshold, quadtree should not subdivide."""
        rng = np.random.default_rng(42)
        image_px = 1000
        image_um = image_px * PIXEL_SIZE

        x = rng.uniform(0, image_um, 1000)
        y = rng.uniform(0, image_um, 1000)

        patches, _, _, _ = compute_density_quadtree_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            pixel_size_um=PIXEL_SIZE,
            x_coords_um=x,
            y_coords_um=y,
            max_transcripts_per_patch=10000,
        )

        # ceil(212.5/100)=3 -> 3x3 = 9 initial patches, no subdivision
        assert len(patches) == 9

    def test_quadtree_dense_region(self):
        """Put 90% of transcripts in one corner, verify subdivision produces more patches."""
        rng = np.random.default_rng(42)
        image_px = 1000
        image_um = image_px * PIXEL_SIZE

        x_sparse = rng.uniform(0, image_um, 100)
        y_sparse = rng.uniform(0, image_um, 100)
        x_dense = rng.uniform(0, image_um * 0.2, 5000)
        y_dense = rng.uniform(0, image_um * 0.2, 5000)
        x = np.concatenate([x_sparse, x_dense])
        y = np.concatenate([y_sparse, y_dense])

        patches, _, _, _ = compute_density_quadtree_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            pixel_size_um=PIXEL_SIZE,
            x_coords_um=x,
            y_coords_um=y,
            max_transcripts_per_patch=500,
            min_tile_width_um=10.0,
            max_depth=4,
        )

        # Should have subdivided beyond the initial 9
        assert len(patches) > 9

    def test_quadtree_max_depth(self):
        """Verify subdivision stops at max_depth: deeper depth -> more patches."""
        rng = np.random.default_rng(42)
        image_px = 1000
        image_um = image_px * PIXEL_SIZE

        x = rng.normal(image_um / 2, 5.0, 10000)
        y = rng.normal(image_um / 2, 5.0, 10000)

        patches_d1, _, _, _ = compute_density_quadtree_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            pixel_size_um=PIXEL_SIZE,
            x_coords_um=x,
            y_coords_um=y,
            max_transcripts_per_patch=10,
            min_tile_width_um=1.0,
            max_depth=1,
        )

        patches_d4, _, _, _ = compute_density_quadtree_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            pixel_size_um=PIXEL_SIZE,
            x_coords_um=x,
            y_coords_um=y,
            max_transcripts_per_patch=10,
            min_tile_width_um=1.0,
            max_depth=4,
        )

        assert len(patches_d4) > len(patches_d1)

    def test_quadtree_min_tile_width(self):
        """Verify subdivision stops at min_tile_width: all cores >= min width."""
        rng = np.random.default_rng(42)
        image_px = 1000
        image_um = image_px * PIXEL_SIZE
        min_tile_um = 30.0

        x = rng.normal(image_um / 2, 5.0, 10000)
        y = rng.normal(image_um / 2, 5.0, 10000)

        patches, _, _, _ = compute_density_quadtree_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            pixel_size_um=PIXEL_SIZE,
            x_coords_um=x,
            y_coords_um=y,
            max_transcripts_per_patch=10,
            min_tile_width_um=min_tile_um,
            max_depth=10,
        )

        for p in patches:
            # Allow 1um rounding tolerance from pixel conversion
            assert p.core_bounds_um.width >= min_tile_um - 1.0, (
                f"Patch {p.patch_id} width {p.core_bounds_um.width:.1f} < min {min_tile_um}"
            )
            assert p.core_bounds_um.height >= min_tile_um - 1.0, (
                f"Patch {p.patch_id} height {p.core_bounds_um.height:.1f} < min {min_tile_um}"
            )


# ---------------------------------------------------------------------------
# Division tests
# ---------------------------------------------------------------------------


class TestDivideTranscripts:
    def test_divide_transcripts_basic(
        self, synthetic_transcripts: Path, tmp_path: Path
    ):
        """Divide synthetic parquet, verify per-patch files are written."""
        output_dir = tmp_path / "output"

        divide_transcripts(
            transcripts_path=synthetic_transcripts,
            output_dir=output_dir,
            image_width_px=6000,
            image_height_px=6000,
            tile_width_um=1000.0,
            overlap_um=50.0,
            balanced=False,
            pixel_size_um=PIXEL_SIZE,
            max_workers=1,
        )

        # Grid metadata should exist
        grid_json = output_dir / "patch_grid.json"
        assert grid_json.exists()

        with open(grid_json) as f:
            metadata = json.load(f)

        # Each patch should have a transcripts.parquet file
        for patch in metadata["patches"]:
            patch_parquet = output_dir / patch["patch_id"] / "transcripts.parquet"
            assert patch_parquet.exists(), f"Missing parquet for {patch['patch_id']}"

    def test_divide_transcripts_coordinates_offset(
        self, synthetic_transcripts: Path, tmp_path: Path
    ):
        """Verify coordinates are offset to patch-local space."""
        output_dir = tmp_path / "output"

        divide_transcripts(
            transcripts_path=synthetic_transcripts,
            output_dir=output_dir,
            image_width_px=6000,
            image_height_px=6000,
            tile_width_um=1000.0,
            overlap_um=50.0,
            balanced=False,
            pixel_size_um=PIXEL_SIZE,
            max_workers=1,
        )

        with open(output_dir / "patch_grid.json") as f:
            metadata = json.load(f)

        for patch_meta in metadata["patches"]:
            patch_parquet = output_dir / patch_meta["patch_id"] / "transcripts.parquet"
            if not patch_parquet.exists():
                continue
            tbl = pq.read_table(str(patch_parquet))
            if tbl.num_rows == 0:
                continue

            gb = patch_meta["global_bounds_um"]
            patch_width = gb["x_max"] - gb["x_min"]
            patch_height = gb["y_max"] - gb["y_min"]

            x_arr = tbl.column("x_location").to_numpy()
            y_arr = tbl.column("y_location").to_numpy()

            # Local coords should be in [0, patch_width) x [0, patch_height)
            assert float(np.min(x_arr)) >= -0.01, (
                f"Patch {patch_meta['patch_id']}: x_min={np.min(x_arr)} < 0"
            )
            assert float(np.max(x_arr)) < patch_width + 0.01, (
                f"Patch {patch_meta['patch_id']}: x_max={np.max(x_arr)} >= {patch_width}"
            )
            assert float(np.min(y_arr)) >= -0.01
            assert float(np.max(y_arr)) < patch_height + 0.01

    def test_divide_transcripts_no_transcript_loss(
        self, synthetic_transcripts: Path, tmp_path: Path
    ):
        """Verify all transcripts appear in at least one patch."""
        output_dir = tmp_path / "output"

        original = pq.read_table(str(synthetic_transcripts))
        original_ids = set(original.column("transcript_id").to_pylist())

        divide_transcripts(
            transcripts_path=synthetic_transcripts,
            output_dir=output_dir,
            image_width_px=6000,
            image_height_px=6000,
            tile_width_um=1000.0,
            overlap_um=50.0,
            balanced=False,
            pixel_size_um=PIXEL_SIZE,
            max_workers=1,
        )

        with open(output_dir / "patch_grid.json") as f:
            metadata = json.load(f)

        found_ids: set[str] = set()
        for patch_meta in metadata["patches"]:
            patch_parquet = output_dir / patch_meta["patch_id"] / "transcripts.parquet"
            if patch_parquet.exists():
                tbl = pq.read_table(str(patch_parquet))
                found_ids.update(tbl.column("transcript_id").to_pylist())

        # Every original transcript must appear in at least one patch
        missing = original_ids - found_ids
        assert len(missing) == 0, f"{len(missing)} transcripts lost during division"


# ---------------------------------------------------------------------------
# Grid metadata JSON roundtrip
# ---------------------------------------------------------------------------


class TestGridMetadataJSON:
    def test_grid_metadata_json_roundtrip(self, tmp_path: Path):
        """Save + load grid metadata preserves all fields."""
        image_px = 1000
        extent = Bounds(0.0, image_px * PIXEL_SIZE, 0.0, image_px * PIXEL_SIZE)
        patches, rows, cols, overlap_px = compute_tilewidth_uniform_grid(
            image_height_px=image_px,
            image_width_px=image_px,
            tile_width_um=100.0,
            overlap_um=50.0,
            pixel_size_um=PIXEL_SIZE,
            transcript_extent_um=extent,
        )

        path = tmp_path / "patch_grid.json"
        save_grid_metadata(
            patches=patches,
            image_height_px=image_px,
            image_width_px=image_px,
            pixel_size_um=PIXEL_SIZE,
            transcript_extent_um=extent,
            grid_rows=rows,
            grid_cols=cols,
            overlap_um=50.0,
            overlap_px=overlap_px,
            grid_type="uniform",
            output_path=path,
        )

        with open(path) as f:
            data = json.load(f)

        assert data["version"] == "1.0"
        assert data["grid_rows"] == rows
        assert data["grid_cols"] == cols
        assert data["overlap_um"] == 50.0
        assert data["overlap_px"] == overlap_px
        assert data["grid_type"] == "uniform"
        assert len(data["patches"]) == len(patches)

        for orig, loaded in zip(patches, data["patches"]):
            assert loaded["patch_id"] == orig.patch_id
            assert loaded["row"] == orig.row
            assert loaded["col"] == orig.col
            assert loaded["global_bounds_px"]["x_min"] == pytest.approx(
                orig.global_bounds_px.x_min
            )
            assert loaded["core_bounds_um"]["y_max"] == pytest.approx(
                orig.core_bounds_um.y_max
            )


# ---------------------------------------------------------------------------
# Merge sparse tiles tests
# ---------------------------------------------------------------------------


def _make_2x2_grid(pixel_size: float = PIXEL_SIZE) -> tuple[list[PatchInfo], int, int]:
    """Build a 2x2 uniform grid on a 1000x1000 pixel image.

    Returns:
        Tuple of (patches, image_width_px, image_height_px).
    """
    image_px = 1000
    overlap_px = int(math.ceil(50.0 / pixel_size))
    patches = _compute_uniform_grid(
        image_height_px=image_px,
        image_width_px=image_px,
        grid_rows=2,
        grid_cols=2,
        overlap_px=overlap_px,
        pixel_size_um=pixel_size,
    )
    return patches, image_px, overlap_px


class TestMergeSparseTiles:
    def test_no_merge_above_threshold(self):
        """All tiles above threshold -- no merging happens."""
        patches, image_px, overlap_px = _make_2x2_grid()
        image_um = image_px * PIXEL_SIZE

        rng = np.random.default_rng(42)
        n = 4000
        x = rng.uniform(0, image_um, n).astype(np.float64)
        y = rng.uniform(0, image_um, n).astype(np.float64)

        merged, merge_count = merge_sparse_tiles(
            patches=patches,
            x_coords_um=x,
            y_coords_um=y,
            overlap_px=overlap_px,
            pixel_size_um=PIXEL_SIZE,
            image_width_px=image_px,
            image_height_px=image_px,
            min_transcripts=100,
        )

        assert merge_count == 0
        assert len(merged) == len(patches)
        merged_ids = {p.patch_id for p in merged}
        original_ids = {p.patch_id for p in patches}
        assert merged_ids == original_ids

    def test_merge_sparse_edge_tile(self):
        """One corner tile has very few transcripts -- it gets merged into a neighbor."""
        patches, image_px, overlap_px = _make_2x2_grid()

        # Put 500 transcripts in each of 3 tiles, 5 in the last tile (row0_col0)
        rng = np.random.default_rng(7)
        # Find the core bounds of each patch to place transcripts correctly
        patch_map = {p.patch_id: p for p in patches}

        sparse_id = patches[0].patch_id  # first tile gets very few transcripts
        xs, ys = [], []
        for p in patches:
            cb = p.core_bounds_um
            n = 5 if p.patch_id == sparse_id else 500
            xs.append(rng.uniform(cb.x_min + 0.1, cb.x_max - 0.1, n))
            ys.append(rng.uniform(cb.y_min + 0.1, cb.y_max - 0.1, n))

        x = np.concatenate(xs)
        y = np.concatenate(ys)

        merged, merge_count = merge_sparse_tiles(
            patches=patches,
            x_coords_um=x,
            y_coords_um=y,
            overlap_px=overlap_px,
            pixel_size_um=PIXEL_SIZE,
            image_width_px=image_px,
            image_height_px=image_px,
            min_transcripts=100,
        )

        assert merge_count == 1
        assert len(merged) == 3

        # The sparse tile should no longer exist as a patch
        merged_ids = {p.patch_id for p in merged}
        assert sparse_id not in merged_ids

        # The absorbing neighbor's bounds should cover the sparse tile's area
        sparse_global = patch_map[sparse_id].global_bounds_um
        absorber = [p for p in merged if p.patch_id != sparse_id]
        # At least one neighbor should now have bounds covering the sparse tile's origin
        covers_sparse = any(
            p.global_bounds_um.x_min <= sparse_global.x_min + 0.01
            and p.global_bounds_um.y_min <= sparse_global.y_min + 0.01
            for p in absorber
        )
        assert covers_sparse, "No merged tile covers the sparse tile's region"

    def test_merge_preserves_all_transcripts(self, tmp_path: Path):
        """After merging, divide_transcripts with merged grid loses no transcripts."""
        image_px = 1000
        image_um = image_px * PIXEL_SIZE

        # Create transcripts: sparse in one corner, dense elsewhere
        rng = np.random.default_rng(33)
        n_sparse = 10
        n_dense = 990

        x_sparse = rng.uniform(0, image_um * 0.1, n_sparse).astype(np.float32)
        y_sparse = rng.uniform(0, image_um * 0.1, n_sparse).astype(np.float32)
        x_dense = rng.uniform(image_um * 0.3, image_um, n_dense).astype(np.float32)
        y_dense = rng.uniform(image_um * 0.3, image_um, n_dense).astype(np.float32)

        n = n_sparse + n_dense
        table = pa.table(
            {
                "transcript_id": pa.array(
                    [f"tx_{i}" for i in range(n)], type=pa.string()
                ),
                "cell_id": pa.array(["UNASSIGNED"] * n, type=pa.string()),
                "overlaps_nucleus": pa.array([0] * n, type=pa.int32()),
                "feature_name": pa.array(
                    [f"gene_{i % 50}" for i in range(n)], type=pa.string()
                ),
                "x_location": pa.array(
                    np.concatenate([x_sparse, x_dense]), type=pa.float32()
                ),
                "y_location": pa.array(
                    np.concatenate([y_sparse, y_dense]), type=pa.float32()
                ),
                "z_location": pa.array(
                    rng.uniform(0, 10, n).astype(np.float32), type=pa.float32()
                ),
                "qv": pa.array(
                    rng.uniform(20, 40, n).astype(np.float32), type=pa.float32()
                ),
            }
        )
        parquet_path = tmp_path / "transcripts.parquet"
        pq.write_table(table, str(parquet_path))

        original_ids = set(table.column("transcript_id").to_pylist())
        output_dir = tmp_path / "output"

        divide_transcripts(
            transcripts_path=parquet_path,
            output_dir=output_dir,
            image_width_px=image_px,
            image_height_px=image_px,
            tile_width_um=100.0,
            overlap_um=10.0,
            balanced=False,
            pixel_size_um=PIXEL_SIZE,
            max_workers=1,
            min_transcripts=50,
        )

        with open(output_dir / "patch_grid.json") as f:
            metadata = json.load(f)

        found_ids: set[str] = set()
        for patch_meta in metadata["patches"]:
            patch_parquet = output_dir / patch_meta["patch_id"] / "transcripts.parquet"
            if patch_parquet.exists():
                tbl = pq.read_table(str(patch_parquet))
                found_ids.update(tbl.column("transcript_id").to_pylist())

        missing = original_ids - found_ids
        assert len(missing) == 0, (
            f"{len(missing)} transcripts lost after merge + divide"
        )

    def test_merge_disabled_with_zero_threshold(self):
        """min_transcripts=0 disables merging regardless of transcript counts."""
        patches, image_px, overlap_px = _make_2x2_grid()

        # Put only 1 transcript per tile -- still no merge with threshold=0
        rng = np.random.default_rng(99)
        xs, ys = [], []
        for p in patches:
            cb = p.core_bounds_um
            xs.append(rng.uniform(cb.x_min + 0.1, cb.x_max - 0.1, 1))
            ys.append(rng.uniform(cb.y_min + 0.1, cb.y_max - 0.1, 1))

        x = np.concatenate(xs)
        y = np.concatenate(ys)

        merged, merge_count = merge_sparse_tiles(
            patches=patches,
            x_coords_um=x,
            y_coords_um=y,
            overlap_px=overlap_px,
            pixel_size_um=PIXEL_SIZE,
            image_width_px=image_px,
            image_height_px=image_px,
            min_transcripts=0,
        )

        assert merge_count == 0
        assert len(merged) == len(patches)

    def test_count_transcripts_per_tile(self):
        """Unit test for _count_transcripts_per_tile with known placement."""
        patches, image_px, _ = _make_2x2_grid()

        # Place 10 transcripts in each patch's core
        rng = np.random.default_rng(11)
        xs, ys = [], []
        expected_per_patch: dict[str, int] = {}
        counts_list = [10, 20, 30, 40]
        for p, n in zip(patches, counts_list):
            cb = p.core_bounds_um
            xs.append(rng.uniform(cb.x_min + 0.1, cb.x_max - 0.1, n))
            ys.append(rng.uniform(cb.y_min + 0.1, cb.y_max - 0.1, n))
            expected_per_patch[p.patch_id] = n

        x = np.concatenate(xs)
        y = np.concatenate(ys)

        counts = _count_transcripts_per_tile(patches, x, y)

        for pid, expected in expected_per_patch.items():
            assert counts[pid] == expected, (
                f"Patch {pid}: expected {expected}, got {counts[pid]}"
            )

    def test_find_adjacent_patches(self):
        """Each tile in a 2x2 grid has exactly 2 neighbors."""
        patches, _, _ = _make_2x2_grid()
        adjacency = _find_adjacent_patches(patches)

        # 2x2 grid: each corner tile touches 2 others (horizontal + vertical)
        for p in patches:
            neighbors = adjacency[p.patch_id]
            assert len(neighbors) == 2, (
                f"Patch {p.patch_id} has {len(neighbors)} neighbors, expected 2: {neighbors}"
            )
