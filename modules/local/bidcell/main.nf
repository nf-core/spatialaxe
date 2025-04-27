process BIDCELL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // container build details
    // Singularity: https://wave.seqera.io/view/builds/bd-a640305ce8418a41_1
    // Docker: https://wave.seqera.io/view/builds/bd-2c6f091afdcdb4c5_1
    container 'docker.io/dongzehe/bidcell:latest'

    input:
    tuple val(meta),  path(transcripts), path(dapi), path(ref), path(pos_markers), path(neg_markers)

    output:
    tuple val(meta), path("${prefix}")                    , emit: outdir
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Output meta needs to correspond to the input used
    prefix = task.ext.prefix ?: "${meta.id}"
    config = task.ext.bidcell_xenium_config ?: ''

    template "bidcell.py"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    cd ${prefix}
    touch affine.csv
    touch dapi_resized.tif
    touch nuclei.tif
    touch all_gene_names.txt
    touch nuclei_cell_type.h5
    mkdir expr_maps
    mkdir cell_gene_matrices
    mkdir model_outputs
    mkdir transcripts_processed.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        bidcell: \$(python3 -c 'import bidcell; print(bidcell.__version__)')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS

    """
}
