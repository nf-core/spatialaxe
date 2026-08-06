process IMAGE_QC_ANALYSIS {
    tag "${meta.id}"
    label 'process_gpu_qc'

    conda "${moduleDir}/environment.yml"
    // Built from environment.yml (see the pipeline Dockerfile assets). Hosted on the
    // author's quay.io namespace for now; to be migrated to the nf-core org before release.
    container "quay.io/dongzehe/image_qc:1.0.0"

    input:
    tuple val(meta), val(parameters), path(input_files)
    path(roi_thresholds_yaml)

    output:
    tuple val(meta), path(outdir), emit: outdir
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('scikit-image'), eval("python3 -c 'import skimage; print(skimage.__version__)'"), topic: versions, emit: versions_skimage

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
    args << "--sample-id '${meta.id}'"

    // ROI thresholds YAML (staged file)
    if (roi_thresholds_yaml.name != 'NO_FILE') {
        args << "--roi-thresholds-yaml '${roi_thresholds_yaml}'"
    }
    // Parse optional parameters from the parameters list
    def param_map = [:]
    if (parameters) {
        parameters.collate(2).each { key, value ->
            if (value != null && value != '') {
                param_map[key] = value
            }
        }
    }

    // Pass stain-names as a single semicolon-separated string
    if (param_map.containsKey('STAIN_NAMES') && param_map['STAIN_NAMES']) {
        def stains_str = param_map['STAIN_NAMES'].toString()
        args << "--stain-names '${stains_str}'"
    }
    else {
        // Module-level default stain/channel names when none are provided
        args << "--stain-names 'DAPI;Boundary (ATP1A1/E-Cadherin/CD45);Interior - RNA (18S);Protein (alphaSMA/Vimentin)'"
    }

    // ROI size parameter
    if (param_map.containsKey('ROI_SIZE') && param_map['ROI_SIZE']) {
        args << "--roi-size ${param_map['ROI_SIZE']}"
    }
    else {
        args << "--roi-size 35"
    }

    // Analysis-tuning flags (--legacy-focus, --no-snr, --snr-*, --save-dapi-maps-tiff,
    // --lap-sigma) come from conf/modules.config ext.args, keeping this module
    // parameter-agnostic.
    """
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    image_qc.py \\
        ${args.join(' \\\n        ')} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outdir = prefix
    """
    mkdir -p "${outdir}/figures"
    touch "${outdir}/image_qc_metrics.json"
    touch "${outdir}/image_qc_metrics.csv"
    touch "${outdir}/versions.yml"
    """
}
