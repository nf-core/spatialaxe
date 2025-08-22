process BAYSOR_CREATE_DATASET {
    tag "$meta.id"
    label 'process_medium'

    container "khersameesh24/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts)
    val(sample_fraction)

    output:
    tuple val(meta), path("sampled_transcripts.csv"), emit: sampled_transcripts
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'create_dataset.py'

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    """
    touch sampled_transcripts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor: 0.7.1
    END_VERSIONS
    """
}
