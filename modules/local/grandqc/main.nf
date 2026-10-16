process GRANDQC {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(image, name: 'source_image.tif'), path(tissue_model, name: 'tissue_model.pth'),
        path(qc_model, name: 'qc_model.pth'), path(script_files),
        val(mpp), val(mpp_model)

    output:
    tuple val(meta), path('grandqc_results'), emit: results
    tuple val("${task.process}"), val('python'),
        eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def mpp_arg = mpp == null ? '' : "--mpp ${mpp}"
    """
    mkdir -p slides grandqc_results
    ln -s "\$(realpath -- source_image.tif)" slides/slide.tif

    python wsi_tis_detect.py \\
        --slide_folder slides \\
        --output_dir grandqc_results \\
        --device cpu \\
        --model_path tissue_model.pth \\
        ${mpp_arg}

    python main.py \\
        --slide_folder slides \\
        --output_dir grandqc_results \\
        --mpp_model ${mpp_model} \\
        --create_geojson Y \\
        --device cpu \\
        --model_path qc_model.pth \\
        ${mpp_arg}
    """

    stub:
    """
    mkdir -p grandqc_results
    touch grandqc_results/GRANDQC_STUB.txt
    """
}
