library(Matrix)
library(rhdf5)

source <- commandArgs(trailingOnly = TRUE)[1]
dir.create("fixture")
path <- "fixture/counts.h5"
if (source == "synthetic") {
    counts <- Matrix(matrix(c(8, 8, 99, 88, 4, 12, 77, 66, 0, 0, 55, 44), nrow = 4), sparse = TRUE)
    genes <- c("gene_1", "gene_2", "protein_1", "control_1")
    cells <- c("mixed", "pure", "empty")
    h5createFile(path)
    h5createGroup(path, "matrix")
    h5createGroup(path, "matrix/features")
    h5write(as.integer(dim(counts)), path, "matrix/shape")
    h5write(counts@i, path, "matrix/indices")
    h5write(counts@p, path, "matrix/indptr")
    h5write(counts@x, path, "matrix/data")
    h5write(cells, path, "matrix/barcodes")
    h5write(genes, path, "matrix/features/id")
    h5write(c("Gene-One", "Gene-Two", "Protein", "Control"), path, "matrix/features/name")
    h5write(c("Gene Expression", "Gene Expression", "Protein Expression", "Negative Control Probe"), path, "matrix/features/feature_type")
    genes <- genes[1:2]
} else {
    stopifnot(file.copy(source, path))
    types <- as.vector(h5read(path, "matrix/features/feature_type"))
    genes <- as.vector(h5read(path, "matrix/features/id"))[types == "Gene Expression"]
    cells <- as.vector(h5read(path, "matrix/barcodes"))
}

reference <- matrix(1, nrow = 2, ncol = length(genes), dimnames = list(c("synthetic_A", "synthetic_B"), genes))
reference[1, seq.int(1, length(genes), by = 2)] <- 3
reference[2, seq.int(2, length(genes), by = 2)] <- 3
weights <- matrix(rep(c(0.75, 0.25), length(cells)), ncol = 2, byrow = TRUE,
    dimnames = list(cells, rownames(reference)))
if (source == "synthetic") weights["pure", ] <- c(0, 1)
write_named <- function(x, key, path) {
    table <- data.frame(rownames(x), x, check.names = FALSE)
    names(table)[1] <- key
    write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_named(reference, "cell_type", "fixture/reference.tsv")
write_named(weights, "cell_id", "fixture/weights.tsv")
write_named(reference[2:1, rev(genes), drop = FALSE], "cell_type", "fixture/reference_shuffled.tsv")
write_named(weights[rev(cells), 2:1, drop = FALSE], "cell_id", "fixture/weights_shuffled.tsv")
write_named(weights[-1, , drop = FALSE], "cell_id", "fixture/weights_missing.tsv")
negative <- weights
negative[1, ] <- c(-0.25, 1.25)
write_named(negative, "cell_id", "fixture/weights_negative.tsv")
ties <- weights
ties[1, ] <- c(0.5, 0.5)
write_named(ties, "cell_id", "fixture/weights_tied.tsv")
write_named(weights * seq_along(cells), "cell_id", "fixture/weights_scaled.tsv")
if (source == "synthetic") {
    symbols <- reference
    colnames(symbols) <- c("Gene-One", "Gene-Two")
    write_named(symbols, "cell_type", "fixture/reference_symbols.tsv")
}
if (source == "synthetic") {
    stopifnot(file.copy(path, "fixture/counts_empty.h5"))
    h5write(rep(0, length(counts@x)), "fixture/counts_empty.h5", "matrix/data")
}
