process SCS_SEGMENT {
    tag "${meta.id}"
    label 'process_high'

    container 'katwre/scs-segment:latest'

    input:
    tuple val(meta), path(scs_input_bgi_tsv), path(morph2d_tif)

    output:
    tuple val(meta), path("spot_prediction_*.txt"), emit: spot_prediction
    tuple val(meta), path("spot2cell_*.txt"), emit: spot2cell
    tuple val(meta), path("spot2nucl_*.txt"), emit: spot2nucl
    tuple val(meta), path("cell_stats_*.txt"), emit: cell_stats
    tuple val(meta), path("cell_masks_*.png"), emit: cell_masks
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SCS_SEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def epochs = task.ext.epochs ?: 100
    def n_neighbor = task.ext.n_neighbor ?: 50
    def r_estimate = task.ext.r_estimate ?: 15
    def bin_size = task.ext.bin_size ?: 3
    def val_ratio = task.ext.val_ratio ?: 0.0625
    def prealigned = task.ext.prealigned != null ? task.ext.prealigned : false
    def align = task.ext.align != null ? task.ext.align : "None"
    def patch_size = task.ext.patch_size ?: 0
    def prealigned_py = prealigned ? "True" : "False"
    def align_py = align == "None" ? "None" : "'${align.toString().replace("'", "\\'")}'"

    """
    mkdir -p data results .matplotlib .numba_cache .cache

    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"

    cat > run_scs.py <<'PY'
import os
from src import scs

scs.segment_cells(
    '${scs_input_bgi_tsv}',
    '${morph2d_tif}',
    prealigned=${prealigned_py},
    align=${align_py},
    patch_size=${patch_size},
    bin_size=${bin_size},
    epochs=${epochs},
    n_neighbor=${n_neighbor},
    r_estimate=${r_estimate},
    val_ratio=${val_ratio},
)
PY

    conda run -n SCS env PYTHONPATH=/app/SCS python run_scs.py

    cp -v results/spot_prediction_*.txt .
    cp -v results/spot2cell_*.txt .
    cp -v results/spot2nucl_*.txt .
    cp -v results/cell_stats_*.txt .
    cp -v results/cell_masks_*.png .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scs: "0.0.0"
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SCS_SEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    """
    touch spot_prediction_0:0:0:0.txt
    touch spot2cell_0:0:0:0.txt
    touch spot2nucl_0:0:0:0.txt
    touch cell_stats_0:0:0:0.txt
    touch cell_masks_0:0:0:0.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scs: "0.0.0"
    END_VERSIONS
    """
}
