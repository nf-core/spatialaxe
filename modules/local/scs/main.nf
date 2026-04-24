process SCS_SEGMENT {
    tag "${meta.id}"
    label 'process_high'

    container 'ghcr.io/katwre/scs:latest'

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

    def epochs = task.ext.epochs ?: 2
    def n_neighbor = task.ext.n_neighbor ?: 10
    def r_estimate = task.ext.r_estimate ?: 12
    def bin_size = task.ext.bin_size ?: 1

    """
    mkdir -p data results

    conda run -n SCS env PYTHONPATH=/app/SCS python - <<'PY'
import os
from src import scs

scs.segment_cells(
    '${scs_input_bgi_tsv}',
    '${morph2d_tif}',
    prealigned=True,
    align=None,
    patch_size=0,
    bin_size=${bin_size},
    epochs=${epochs},
    n_neighbor=${n_neighbor},
    r_estimate=${r_estimate},
)
PY

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
