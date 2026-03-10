process SEGGER_CREATE_DATASET {
    tag "${meta.id}"
    label 'process_high'
    maxForks params.restrict_concurrency ? 1 : 0

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(base_dir)

    output:
    tuple val(meta), path("${prefix}/"), emit: datasetdir
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    def script_path = "/workspace/segger_dev/src/segger/cli/create_dataset_fast.py"
    prefix = task.ext.prefix ?: "${meta.id}"

    // check for platform values
    if (!(params.format in ['xenium'])) {
        error("${params.format} is an invalid platform type.")
    }

    """
    # Set numba cache directory to avoid caching issues in container
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    # Create local bundle directory with symlinks to all original files
    # This is necessary because input files from S3/Fusion are read-only
    # Use absolute paths to avoid broken relative symlinks
    mkdir -p bundle_local
    for item in ${base_dir}/*; do
        # Resolve to absolute path (follow any symlinks)
        abs_path=\$(readlink -f "\$item" 2>/dev/null || realpath "\$item" 2>/dev/null || echo "\$item")
        basename=\$(basename "\$item")
        ln -sf "\$abs_path" "bundle_local/\$basename"
    done

    # Segger expects nucleus_boundaries.parquet but Xenium bundles have cell_boundaries.parquet
    # Create the symlink if nucleus_boundaries doesn't exist but cell_boundaries does
    if [ ! -e "bundle_local/nucleus_boundaries.parquet" ] && [ -e "bundle_local/cell_boundaries.parquet" ]; then
        echo "Creating nucleus_boundaries.parquet symlink from cell_boundaries.parquet"
        cell_bounds_path=\$(readlink -f "bundle_local/cell_boundaries.parquet" 2>/dev/null || realpath "bundle_local/cell_boundaries.parquet" 2>/dev/null)
        ln -sf "\$cell_bounds_path" bundle_local/nucleus_boundaries.parquet
    fi

    # List bundle contents for debugging
    echo "Bundle contents:"
    ls -la bundle_local/

    # Fix: Add parquet column statistics for segger
    echo "Adding statistics to parquet files..."
    python3 - << 'PYEOF'
import pyarrow.parquet as pq
import os

def add_stats(inp, out):
    if not os.path.exists(inp):
        print(f"  Skip {inp}")
        return
    t = pq.read_table(inp)
    pq.write_table(t, out, write_statistics=True, compression='snappy')
    print(f"  Done {os.path.basename(inp)} ({len(t)} rows)")

os.makedirs('bundle_stats', exist_ok=True)
for f in ['transcripts.parquet', 'nucleus_boundaries.parquet']:
    add_stats(f'bundle_local/{f}', f'bundle_stats/{f}')

for item in os.listdir('bundle_local'):
    s, d = f'bundle_local/{item}', f'bundle_stats/{item}'
    if not os.path.exists(d):
        os.symlink(os.path.realpath(s), d)
print("Done")

# Debug: Check overlaps_nucleus column data
print("")
print("=== Debugging overlaps_nucleus data ===")
import pyarrow.compute as pc

tx = pq.read_table('bundle_stats/transcripts.parquet')
bd = pq.read_table('bundle_stats/nucleus_boundaries.parquet')

if 'overlaps_nucleus' in tx.column_names:
    col = tx.column('overlaps_nucleus')
    print(f"overlaps_nucleus dtype: {col.type}")
    unique_vals = pc.unique(col)
    print(f"overlaps_nucleus unique values: {unique_vals.to_pylist()[:10]}")
    val_counts = pc.value_counts(col)
    print(f"overlaps_nucleus value_counts: {val_counts.to_pylist()}")
else:
    print("WARNING: overlaps_nucleus column NOT FOUND in transcripts.parquet")

# Check cell_id overlap between transcripts and boundaries
if 'cell_id' in tx.column_names and 'cell_id' in bd.column_names:
    tx_cells = set(pc.unique(tx.column('cell_id')).to_pylist())
    bd_cells = set(pc.unique(bd.column('cell_id')).to_pylist())
    overlap = tx_cells & bd_cells
    print("")
    print(f"Transcripts unique cell_ids: {len(tx_cells)}")
    print(f"Boundaries unique cell_ids: {len(bd_cells)}")
    print(f"Overlapping cell_ids: {len(overlap)}")

print("=== End Debug ===")
PYEOF
    ls -la bundle_stats/

    python3 ${script_path} \\
        --base_dir bundle_stats \\
        --data_dir ${prefix} \\
        --sample_type ${params.format} \\
        --tile_width ${params.tile_width} \\
        --tile_height ${params.tile_height} \\
        --n_workers ${task.cpus} \\
        ${args}

    # Verify tiles were created and show distribution
    echo "Dataset split (before fix):"
    echo "  train_tiles: \$(ls ${prefix}/train_tiles/processed/ 2>/dev/null | wc -l) files"
    echo "  val_tiles: \$(ls ${prefix}/val_tiles/processed/ 2>/dev/null | wc -l) files"
    echo "  test_tiles: \$(ls ${prefix}/test_tiles/processed/ 2>/dev/null | wc -l) files"

    # Workaround: segger commit 0787167 has a bug where all tiles go to test_tiles
    # regardless of test_prob/val_prob settings. Move ONLY trainable tiles (those with
    # edge_label_index) from test_tiles to train_tiles.
    # Tiles without tx-belongs-bd edges don't have edge_label_index and cannot be used for training.
    train_count=\$(ls ${prefix}/train_tiles/processed/ 2>/dev/null | wc -l)
    test_count=\$(ls ${prefix}/test_tiles/processed/ 2>/dev/null | wc -l)

    if [ "\$train_count" -eq 0 ] && [ "\$test_count" -gt 0 ]; then
        echo "Applying workaround: filtering trainable tiles from test_tiles (segger split bug)"
        export SEGGER_PREFIX="${prefix}"
        python3 - << 'PYEOF'
import torch
import os
import shutil

prefix = os.environ['SEGGER_PREFIX']
test_dir = f"{prefix}/test_tiles/processed"
train_dir = f"{prefix}/train_tiles/processed"

moved = 0
skipped = 0

for f in os.listdir(test_dir):
    if not f.endswith('.pt'):
        continue
    fpath = os.path.join(test_dir, f)
    try:
        tile = torch.load(fpath, weights_only=False)
        edge_store = tile['tx', 'belongs', 'bd']
        # Check if edge_label_index exists and has data
        if hasattr(edge_store, 'edge_label_index') and edge_store.edge_label_index.numel() > 0:
            shutil.move(fpath, os.path.join(train_dir, f))
            moved += 1
        else:
            skipped += 1
    except Exception as e:
        print(f"Warning: Could not process {f}: {e}")
        skipped += 1

print(f"Moved {moved} trainable tiles to train_tiles")
print(f"Skipped {skipped} test-only tiles (no edge_label_index)")
PYEOF
    fi

    echo "Dataset split (after fix):"
    echo "  train_tiles: \$(ls ${prefix}/train_tiles/processed/ 2>/dev/null | wc -l) files"
    echo "  val_tiles: \$(ls ${prefix}/val_tiles/processed/ 2>/dev/null | wc -l) files"
    echo "  test_tiles: \$(ls ${prefix}/test_tiles/processed/ 2>/dev/null | wc -l) files"

    train_tiles_dir="${prefix}/train_tiles/processed"
    if [ ! -d "\$train_tiles_dir" ] || [ -z "\$(ls -A \$train_tiles_dir 2>/dev/null)" ]; then
        echo "ERROR: No trainable tiles were created in \$train_tiles_dir"
        echo "This usually means no transcripts overlap with nucleus boundaries in the dataset."
        echo "Check if the Xenium bundle contains valid overlaps_nucleus data in transcripts.parquet."
        exit 1
    fi
    echo "Successfully created \$(ls \$train_tiles_dir | wc -l) trainable tiles"

    # Workaround: Segger's get_polygon_props() produces NaN boundary features (bd.x)
    # when polygon geometries have zero area or index misalignment during GeoDataFrame
    # construction. Replace NaN bd.x with zeros so BCEWithLogitsLoss doesn't propagate NaN.
    export SEGGER_PREFIX="${prefix}"
    python3 - << 'PYEOF'
import torch
import os

prefix = os.environ['SEGGER_PREFIX']
fixed = 0
total = 0

for split in ['train_tiles', 'test_tiles', 'val_tiles']:
    tile_dir = f"{prefix}/{split}/processed"
    if not os.path.isdir(tile_dir):
        continue
    for f in os.listdir(tile_dir):
        if not f.endswith('.pt'):
            continue
        total += 1
        fpath = os.path.join(tile_dir, f)
        tile = torch.load(fpath, weights_only=False)
        bd_x = tile['bd'].x
        if bd_x.isnan().any():
            tile['bd'].x = torch.nan_to_num(bd_x, nan=0.0)
            torch.save(tile, fpath)
            fixed += 1

print(f"Fixed NaN bd.x in {fixed}/{total} tiles")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch "${prefix}/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """
}
