process QUARTO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Shared report renderer. In this pipeline each QC subworkflow overrides the
    // container in conf/modules.config to render inside its own QC image (both of
    // which include Quarto); this default is one of those images.
    container "quay.io/dongzehe/image_qc:1.0.0"

    input:
    tuple val(meta), path(indir)
    tuple val(meta2), path(qmd_file)
    path(roi_thresholds_yaml)

    output:
    tuple val(meta), path("${html_name}.${format}"), emit: report
    tuple val("${task.process}"), val('quarto'), eval("quarto --version"), topic: versions, emit: versions_quarto

    when:
    task.ext.when == null || task.ext.when

    script:
    // Process-scoped (no `def`): required so output: can resolve ${html_name}
    html_name = task.ext.prefix ?: qmd_file.baseName
    def quarto_args = task.ext.quarto_args ?: '--embed-resources --standalone'
    format = task.ext.quarto_format ?: 'html'

    // Shell-safe single-quoted strings for quarto -P (sample id / paths may contain quotes).
    // Inline the sq escape rather than defining a `def sq = { ... }` closure:
    // Nextflow 26's strict parser can't resolve closure calls inside GString
    // interpolations (`${sq(...)}`), so we pre-compute each escaped value.
    def bundleStr = (meta.containsKey('samplesheet_xenium_bundle') && meta.samplesheet_xenium_bundle != null) ? meta.samplesheet_xenium_bundle.toString() : ''
    // Expected publish dir for this step (matches the QC withName blocks in conf/modules.config:
    // ${params.outdir}/${params.mode}/qc/<prefix>). Do not reference html_name here —
    // Nextflow's parser treats it as redefining html_name vs output:
    def pubStr = bundleStr ? "${params.outdir}/${params.mode}/qc/${task.ext.prefix ?: qmd_file.baseName}" : ''
    def sqSample = meta.id.toString().replace("'", "'\\''")
    def sqBundle = bundleStr.replace("'", "'\\''")
    def sqPub = pubStr.replace("'", "'\\''")
    // Image QC report only: pass staged YAML so Quarto finds it in the container work dir
    def roiYamlLine = ''
    if (task.ext.pass_roi_thresholds_yaml && roi_thresholds_yaml.size() > 0) {
        def sqRoi = roi_thresholds_yaml.name.replace("'", "'\\''")
        roiYamlLine = "        -P ROI_THRESHOLDS_YAML:'${sqRoi}' \\\n"
    }

    """
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    quarto render ${qmd_file} \\
        --to ${format} \\
        ${quarto_args} \\
        -P INDIR:${indir} \\
${roiYamlLine}        -P SAMPLE_NAME:'${sqSample}' \\
        -P XENIUM_BUNDLE:'${sqBundle}' \\
        -P SAMPLE_PUBLISHED_OUTDIR:'${sqPub}' \\
        --output ${html_name}.${format}
    """

    stub:
    html_name = task.ext.prefix ?: qmd_file.baseName
    format = task.ext.quarto_format ?: 'html'
    """
    touch "${html_name}.${format}"
    """
}
