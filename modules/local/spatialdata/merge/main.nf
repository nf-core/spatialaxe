process SPATIALDATA_MERGE {
    tag "$meta.id"
    label 'process_low'

    container "heylf/spatialdata:0.2.6"

    input:
    tuple val(meta), path(ref_bundle, stageAs: "*")
    path(add_bundle, stageAs: "*")

    output:
    tuple val(meta), path("spatialdata_spatialxe"), emit: spatialxe_bundle
    path("versions.yml")                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "SPATIALDATA_WRITE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''

    template 'merge.py'

    stub:
    """
    mkdir -p "spatialdata_spatialxe/"
    touch spatialdata_spatialxe/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata: \$(echo \$( python -c "import spatialdata; print(spatialdata.__version__)" 2>&1) )
    END_VERSIONS
    """

}
