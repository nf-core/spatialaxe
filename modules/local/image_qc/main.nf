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

    // --roi-size comes from conf/modules.config ext.args (params.image_qc_roi_size),
    // keeping this module parameter-agnostic.

    // GPU cap: pass --num-gpus to the script only when the user set params.num_gpus.
    // Mirror modules/local/segger/train pattern -- constrain CUDA_VISIBLE_DEVICES so
    // the process never grabs more GPUs than requested (accelerator directive from
    // conf/base.config process_gpu_qc already requests this many devices).
    def num_gpus = params.num_gpus
    def cuda_visible = ''
    if (num_gpus != null) {
        args << "--num-gpus ${num_gpus}"
        cuda_visible = (num_gpus as int) > 0
            ? "export CUDA_VISIBLE_DEVICES=" + (0..<(num_gpus as int)).join(',')
            : "export CUDA_VISIBLE_DEVICES="
    }

    // Analysis-tuning flags (--legacy-focus, --no-snr, --snr-*, --save-dapi-maps-tiff,
    // --lap-sigma) come from conf/modules.config ext.args, keeping this module
    // parameter-agnostic.
    """
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"
    ${cuda_visible}

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
