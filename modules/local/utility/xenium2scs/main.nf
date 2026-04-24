process XENIUM2SCS {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/katwre/scs-segment:latest"

    input:
    tuple val(meta), path(transcripts_parquet), path(morphology_image)

    output:
    tuple val(meta), path("${prefix}/scs_input_bgi.tsv"), emit: scs_input_bgi_tsv
    tuple val(meta), path("${prefix}/morph2d.tif"), emit: morph2d_tif
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("XENIUM2SCS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    template('xenium2scs.py')

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("XENIUM2SCS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    printf 'geneID\tx\ty\tMIDCounts\nTEST\t0\t0\t1\n' > ${prefix}/scs_input_bgi.tsv

    python - <<'PY'
import numpy as np
import tifffile

img = np.zeros((16, 16), dtype=np.uint16)
tifffile.imwrite('${prefix}/morph2d.tif', img)
PY

    printf '"${task.process}":\n    xenium2scs: "1.0.0"\n' > versions.yml
    """
}
