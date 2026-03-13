process XENIUM2SCS {
    tag "$meta.id"
    label 'process_low'

    container "khersameesh24/spatialdata:0.2.6"

    input:
    tuple val(meta), path(transcripts_parquet), path(morphology_image), path(experiment_xenium)

    output:
    tuple val(meta), path("${prefix}/scs_input_bgi.tsv"), emit: scs_input_bgi_tsv
    tuple val(meta), path("${prefix}/morph2d.tif"), emit: morph2d_tif
    tuple val(meta), path("${prefix}/xenium2scs_metrics.tsv"), emit: metrics
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

    cat <<'EOF' > ${prefix}/scs_input.tsv
    geneID\trow\tcolumn\tcounts
    TEST\t0\t0\t1
    EOF

    cat <<'EOF' > ${prefix}/scs_input_bgi.tsv
    geneID\tx\ty\tMIDCounts
    TEST\t0\t0\t1
    EOF

    python - <<'PY'
import numpy as np
import tifffile

img = np.zeros((16, 16), dtype=np.uint16)
tifffile.imwrite('${prefix}/morph2d.tif', img)
PY

    cat <<'EOF' > ${prefix}/xenium2scs_metrics.tsv
    metric\tvalue
    n_rows\t1
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xenium2scs: "1.0.0"
    END_VERSIONS
    """
}
