process SPLIT_TEST_INPUTS {
    container 'docker.io/library/split:0.3.0'

    input:
    path counts
    path generator

    output:
    path 'fixture', emit: fixture

    script:
    def source = counts ? counts : 'synthetic'
    """
    Rscript ${generator} ${source}
    """
}
