nextflow.enable.dsl = 2

include { GRANDQC } from './modules/local/grandqc/main'

workflow {
    grandqc_input = Channel.of([
        [id: 'ovarian'],
        file(params.image, checkIfExists: true),
        file(params.tissue_model, checkIfExists: true),
        file(params.qc_model, checkIfExists: true),
        files("${projectDir}/updated_scripts/*.py"),
        params.grandqc_mpp,
        params.grandqc_mpp_model ?: 1.5
    ])

    GRANDQC(grandqc_input)
}
