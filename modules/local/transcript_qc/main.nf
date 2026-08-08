process TRANSCRIPT_QC_PROCESSING {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Built from environment.yml (see the pipeline Dockerfile assets). Hosted on the
    // author's quay.io namespace for now; to be migrated to the nf-core org before release.
    container "quay.io/dongzehe/transcript_qc:1.0.0"

    input:
    tuple val(meta), val(parameters), path(input_files)

    output:
    tuple val(meta), path(outdir), emit: outdir
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('scanpy'), eval("python3 -c 'import scanpy; print(scanpy.__version__)'"), topic: versions, emit: versions_scanpy
    tuple val("${task.process}"), val('anndata'), eval("python3 -c 'import anndata; print(anndata.__version__)'"), topic: versions, emit: versions_anndata

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outdir = prefix

    // Convert parameters to script arguments
    def args = []

    // Required parameters
    args << "--xenium-bundle-dir '${input_files[0]}'"
    args << "--outdir '${outdir}'"
    args << "--threads ${task.cpus}"
    args << "--task-process '${task.process}'"

    // Parse optional parameters from the parameters list
    def param_map = [:]
    if (parameters) {
        parameters.collate(2).each { key, value ->
            if (value != null && value != '') {
                param_map[key] = value
            }
        }
    }

    // Add optional parameters
    if (param_map.containsKey('NON_GENE_PREFIX') && param_map['NON_GENE_PREFIX']) {
        def prefixes = param_map['NON_GENE_PREFIX'].toString().split(';').collect { "'${it.trim()}'" }.join(' ')
        args << "--non-gene-prefix ${prefixes}"
    }

    if (param_map.containsKey('STAIN_NAMES') && param_map['STAIN_NAMES']) {
        def stains = param_map['STAIN_NAMES'].toString().split(';').collect { "'${it.trim()}'" }.join(' ')
        args << "--stain-names ${stains}"
    }

    if (param_map.containsKey('NUM_ROW_GROUPS') && param_map['NUM_ROW_GROUPS']) {
        args << "--num-row-groups ${param_map['NUM_ROW_GROUPS']}"
    }

    """
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    transcript_qc_processing.py \\
        ${args.join(' \\\n        ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outdir = prefix
    """
    mkdir -p "${outdir}/figures"
    touch "${outdir}/transcript_qc_metrics.json"
    touch "${outdir}/versions.yml"
    """
}
