process COLLECT_PIPELINE_QC {
    tag "pipeline_qc"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::pyyaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/9409ce399922a5746bea1b7df5668c3d1d79b9af49a15950d9818c4fe45ac749/data' :
        'community.wave.seqera.io/library/pandas_procs_pyarrow:d8f882b65dfea451' }"

    input:
    path(qc_jsons)

    output:
    path("mqc_output/*_mqc.yaml"), emit: mqc_custom_content
    path("all_qc_records.json"),   emit: merged_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p mqc_output

    emit_pipeline_qc.py \\
        --qc-files ${qc_jsons} \\
        --output-dir mqc_output \\
        --merged-output all_qc_records.json
    """

    stub:
    """
    mkdir -p mqc_output
    echo '{"id":"spatialaxe_pipeline_summary","section_name":"Pipeline Summary","data":{}}' > mqc_output/spatialaxe_pipeline_summary_mqc.yaml
    echo '[]' > all_qc_records.json
    """
}
