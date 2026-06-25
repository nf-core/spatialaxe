process PARQUET_TO_CSV {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c8f3a2de4e2aa94500011f7d7d09276e9b6f2d79ee8737c9098fe22d4649bc/data' :
        'community.wave.seqera.io/library/sopa_procps-ng_pyarrow:c9ce8cd2ede79d72' }"

    input:
    tuple val(meta), path(parquet)

    output:
    tuple val(meta), path("transcripts.csv"), emit: csv
    tuple val("${task.process}"), val('pyarrow'), eval("python3 -c 'import pyarrow; print(pyarrow.__version__)'"), topic: versions, emit: versions_pyarrow

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import pyarrow.parquet as pq
import pyarrow.csv as pa_csv
t = pq.read_table('${parquet}')
pa_csv.write_csv(t, 'transcripts.csv')
"
    """

    stub:
    """
    touch transcripts.csv
    """
}
