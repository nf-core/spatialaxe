#!/usr/bin/env python3
"""
Run segger create_dataset with spatialaxe-specific preprocessing and workarounds.

Wraps segger's create_dataset_fast.py with:
  - bundle_local symlink prep (handles read-only S3/Fusion mounts)
  - parquet column statistics (segger needs these)
  - WORKAROUND: filter trainable tiles from test_tiles when segger commit 0787167 mis-splits
  - WORKAROUND: replace NaN bd.x with zeros after get_polygon_props produces NaN

Each WORKAROUND should be removable when the upstream segger bug is fixed.
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# imports for actual work (used in functions below)
import pyarrow.parquet as pq
import pyarrow.compute as pc
import torch


SEGGER_CLI = "/workspace/segger_dev/src/segger/cli/create_dataset_fast.py"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bundle-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--sample-type", required=True, choices=["xenium"])
    p.add_argument("--tile-width", type=int, required=True)
    p.add_argument("--tile-height", type=int, required=True)
    p.add_argument("--n-workers", type=int, required=True)
    # remaining args forwarded to segger CLI
    args, extra = p.parse_known_args()
    return args, extra


def prepare_bundle(bundle_dir):
    """Create local bundle dir with absolute symlinks (S3/Fusion read-only-safe)."""
    Path("bundle_local").mkdir(exist_ok=True)
    for item in Path(bundle_dir).iterdir():
        try:
            abs_path = item.resolve()
        except Exception:
            abs_path = item
        target = Path("bundle_local") / item.name
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(abs_path)

    # Segger expects nucleus_boundaries.parquet but Xenium bundles have cell_boundaries.parquet
    nb = Path("bundle_local/nucleus_boundaries.parquet")
    cb = Path("bundle_local/cell_boundaries.parquet")
    if not nb.exists() and cb.exists():
        print(
            "Creating nucleus_boundaries.parquet symlink from cell_boundaries.parquet"
        )
        nb.symlink_to(cb.resolve())

    print("Bundle contents:")
    for item in sorted(Path("bundle_local").iterdir()):
        print(f"  {item.name}")


def add_parquet_stats():
    """Rewrite key parquet files with column statistics (segger requires them)."""
    Path("bundle_stats").mkdir(exist_ok=True)
    for fname in ["transcripts.parquet", "nucleus_boundaries.parquet"]:
        src = Path("bundle_local") / fname
        dst = Path("bundle_stats") / fname
        if not src.exists():
            print(f"  Skip {src}")
            continue
        t = pq.read_table(str(src))
        pq.write_table(t, str(dst), write_statistics=True, compression="snappy")
        print(f"  Done {fname} ({len(t)} rows)")

    # Symlink everything else from bundle_local into bundle_stats
    for item in Path("bundle_local").iterdir():
        dst = Path("bundle_stats") / item.name
        if not dst.exists():
            dst.symlink_to(item.resolve())

    # Debug: check overlaps_nucleus column in transcripts
    print("\n=== Debugging overlaps_nucleus data ===")
    tx = pq.read_table("bundle_stats/transcripts.parquet")
    bd = pq.read_table("bundle_stats/nucleus_boundaries.parquet")
    if "overlaps_nucleus" in tx.column_names:
        col = tx.column("overlaps_nucleus")
        print(f"overlaps_nucleus dtype: {col.type}")
        unique_vals = pc.unique(col)
        print(f"overlaps_nucleus unique values: {unique_vals.to_pylist()[:10]}")
        val_counts = pc.value_counts(col)
        print(f"overlaps_nucleus value_counts: {val_counts.to_pylist()}")
    else:
        print("WARNING: overlaps_nucleus column NOT FOUND in transcripts.parquet")

    if "cell_id" in tx.column_names and "cell_id" in bd.column_names:
        tx_cells = set(pc.unique(tx.column("cell_id")).to_pylist())
        bd_cells = set(pc.unique(bd.column("cell_id")).to_pylist())
        overlap = tx_cells & bd_cells
        print(f"Transcripts unique cell_ids: {len(tx_cells)}")
        print(f"Boundaries unique cell_ids: {len(bd_cells)}")
        print(f"Overlapping cell_ids: {len(overlap)}")
    print("=== End Debug ===\n")


def run_segger_cli(args, extra):
    cmd = [
        "python3",
        SEGGER_CLI,
        "--base_dir",
        "bundle_stats",
        "--data_dir",
        args.output_dir,
        "--sample_type",
        args.sample_type,
        "--tile_width",
        str(args.tile_width),
        "--tile_height",
        str(args.tile_height),
        "--n_workers",
        str(args.n_workers),
        *extra,
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def filter_trainable_tiles_if_needed(prefix):
    """
    WORKAROUND: segger commit 0787167 has a bug where all tiles end up in test_tiles
    regardless of test_prob/val_prob settings. Move ONLY trainable tiles (those with
    edge_label_index) from test_tiles to train_tiles.

    Remove this function once segger >= 0.1.x is bumped with the upstream fix.
    """
    train_dir = Path(prefix) / "train_tiles" / "processed"
    test_dir = Path(prefix) / "test_tiles" / "processed"
    val_dir = Path(prefix) / "val_tiles" / "processed"

    train_count = len(list(train_dir.iterdir())) if train_dir.exists() else 0
    test_count = len(list(test_dir.iterdir())) if test_dir.exists() else 0
    val_count = len(list(val_dir.iterdir())) if val_dir.exists() else 0
    print(
        f"Dataset split (before fix): train={train_count} val={val_count} test={test_count}"
    )

    if train_count == 0 and test_count > 0:
        print(
            "Applying workaround: filtering trainable tiles from test_tiles (segger split bug)"
        )
        moved = 0
        skipped = 0
        for tile_path in list(test_dir.iterdir()):
            if not tile_path.name.endswith(".pt"):
                continue
            try:
                tile = torch.load(str(tile_path), weights_only=False)
                edge_store = tile["tx", "belongs", "bd"]
                if (
                    hasattr(edge_store, "edge_label_index")
                    and edge_store.edge_label_index.numel() > 0
                ):
                    shutil.move(str(tile_path), str(train_dir / tile_path.name))
                    moved += 1
                else:
                    skipped += 1
            except Exception as e:
                print(f"Warning: Could not process {tile_path.name}: {e}")
                skipped += 1
        print(f"Moved {moved} trainable tiles to train_tiles")
        print(f"Skipped {skipped} test-only tiles (no edge_label_index)")

    train_count = len(list(train_dir.iterdir())) if train_dir.exists() else 0
    test_count = len(list(test_dir.iterdir())) if test_dir.exists() else 0
    val_count = len(list(val_dir.iterdir())) if val_dir.exists() else 0
    print(
        f"Dataset split (after fix): train={train_count} val={val_count} test={test_count}"
    )

    if train_count == 0:
        print(f"ERROR: No trainable tiles were created in {train_dir}", file=sys.stderr)
        print(
            "This usually means no transcripts overlap with nucleus boundaries in the dataset.",
            file=sys.stderr,
        )
        print(
            "Check if the Xenium bundle contains valid overlaps_nucleus data in transcripts.parquet.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"Successfully created {train_count} trainable tiles")


def fix_bd_x_nan(prefix):
    """
    WORKAROUND: segger's get_polygon_props() produces NaN boundary features (bd.x)
    when polygon geometries have zero area or index misalignment during GeoDataFrame
    construction. Replace NaN bd.x with zeros so BCEWithLogitsLoss doesn't propagate NaN.

    Remove this function once segger >= 0.1.x is bumped with the upstream fix.
    """
    fixed = 0
    total = 0
    for split in ["train_tiles", "test_tiles", "val_tiles"]:
        tile_dir = Path(prefix) / split / "processed"
        if not tile_dir.is_dir():
            continue
        for tile_path in tile_dir.iterdir():
            if not tile_path.name.endswith(".pt"):
                continue
            total += 1
            tile = torch.load(str(tile_path), weights_only=False)
            bd_x = tile["bd"].x
            if bd_x.isnan().any():
                tile["bd"].x = torch.nan_to_num(bd_x, nan=0.0)
                torch.save(tile, str(tile_path))
                fixed += 1
    print(f"Fixed NaN bd.x in {fixed}/{total} tiles")


def main():
    args, extra = parse_args()

    # Ensure numba cache dir is writable (env var should be set by caller, but belt-and-suspenders)
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.getcwd(), ".numba_cache"))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

    prepare_bundle(args.bundle_dir)
    print("Adding statistics to parquet files...")
    add_parquet_stats()

    # Sanity-check bundle_stats
    print("bundle_stats contents:")
    for item in sorted(Path("bundle_stats").iterdir()):
        print(f"  {item.name}")

    run_segger_cli(args, extra)

    filter_trainable_tiles_if_needed(args.output_dir)
    fix_bd_x_nan(args.output_dir)


if __name__ == "__main__":
    main()
