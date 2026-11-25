#!/usr/bin/env Rscript

library(Matrix)
library(rhdf5)
library(SPLIT)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || length(args) %% 2 != 0)
  stop("Options must be supplied as --name value pairs.", call. = FALSE)
keys <- args[seq.int(1, length(args), by = 2)]
allowed <- c(
  "--counts",
  "--reference",
  "--weights",
  "--outdir",
  "--gene-key",
  "--chunk-size"
)
if (anyDuplicated(keys) || any(!keys %in% allowed))
  stop("Unknown or repeated option.", call. = FALSE)
options <- setNames(as.list(args[seq.int(2, length(args), by = 2)]), keys)
required <- c("--counts", "--reference", "--weights", "--outdir")
if (!all(required %in% keys))
  stop(
    "Required options: --counts, --reference, --weights, --outdir.",
    call. = FALSE
  )
gene_key <- if (is.null(options[["--gene-key"]])) "id" else
  options[["--gene-key"]]
if (!gene_key %in% c("id", "name"))
  stop("--gene-key must be id or name.", call. = FALSE)
# SPLIT creates dense cell-by-gene intermediates even for sparse count inputs.
chunk_size <- if (is.null(options[["--chunk-size"]])) 1000 else
  suppressWarnings(as.numeric(options[["--chunk-size"]]))
if (
  !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size)
) {
  stop("--chunk-size must be a positive integer.", call. = FALSE)
}

check_ids <- function(ids, description) {
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(description, " must be nonempty and unique.", call. = FALSE)
  }
}
read_table <- function(path, key) {
  table <- read.delim(
    path,
    check.names = FALSE,
    colClasses = "character",
    na.strings = NULL,
    quote = "",
    comment.char = "",
    row.names = NULL
  )
  if (ncol(table) < 2 || names(table)[1] != key) {
    stop(
      path,
      " must start with a ",
      key,
      " column followed by numeric columns.",
      call. = FALSE
    )
  }
  ids <- table[[1]]
  check_ids(ids, paste(path, "row identifiers"))
  check_ids(names(table)[-1], paste(path, "column identifiers"))
  values <- as.matrix(table[-1])
  suppressWarnings(storage.mode(values) <- "double")
  if (any(!is.finite(values)) || any(values < 0)) {
    stop(
      path,
      " must contain finite, nonnegative numeric values.",
      call. = FALSE
    )
  }
  rownames(values) <- ids
  values
}

path <- options[["--counts"]]
feature_types <- as.vector(h5read(path, "matrix/features/feature_type"))
feature_ids <- as.vector(h5read(path, "matrix/features/id"))
feature_names <- as.vector(h5read(path, "matrix/features/name"))
cells <- as.vector(h5read(path, "matrix/barcodes"))
shape <- as.integer(h5read(path, "matrix/shape"))
if (
  length(shape) != 2 ||
    shape[1] != length(feature_types) ||
    shape[1] != length(feature_ids) ||
    shape[1] != length(feature_names) ||
    shape[2] != length(cells)
) {
  stop(
    "HDF5 matrix dimensions do not match the feature and cell identifiers.",
    call. = FALSE
  )
}
check_ids(cells, "Cell identifiers")
counts <- sparseMatrix(
  i = as.integer(h5read(path, "matrix/indices")) + 1L,
  p = as.integer(h5read(path, "matrix/indptr")),
  x = as.numeric(h5read(path, "matrix/data")),
  dims = shape
)
if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
  stop("Counts must be finite and nonnegative.", call. = FALSE)
}
rna <- feature_types == "Gene Expression"
features <- data.frame(
  id = feature_ids[rna],
  name = feature_names[rna],
  type = feature_types[rna]
)
genes <- features[[gene_key]]
check_ids(genes, paste("RNA feature", gene_key, "identifiers"))
counts <- counts[rna, , drop = FALSE]
dimnames(counts) <- list(genes, cells)
if (nrow(counts) < 2 || ncol(counts) < 2) {
  stop(
    "SPLIT 0.3.0 requires at least two RNA features and two cells.",
    call. = FALSE
  )
}
if (nnzero(counts) == 0) {
  stop(
    "SPLIT 0.3.0 cannot purify an entirely zero RNA count matrix.",
    call. = FALSE
  )
}

reference <- read_table(options[["--reference"]], "cell_type")
weights <- read_table(options[["--weights"]], "cell_id")
if (!setequal(cells, rownames(weights))) {
  stop(
    "Weights must contain exactly the cells in the count matrix.",
    call. = FALSE
  )
}
if (!setequal(rownames(reference), colnames(weights))) {
  stop("Reference and weight cell types must match.", call. = FALSE)
}
if (!all(genes %in% colnames(reference))) {
  stop(
    "Reference must cover every RNA feature selected by --gene-key.",
    call. = FALSE
  )
}
if (nrow(reference) < 2)
  stop("SPLIT 0.3.0 requires at least two reference cell types.", call. = FALSE)
reference <- reference[, genes, drop = FALSE]
weights <- weights[cells, rownames(reference), drop = FALSE]
totals <- rowSums(weights)
if (any(!is.finite(totals)) || any(totals <= 0))
  stop("Each cell must have positive total weight.", call. = FALSE)
weights <- weights / totals
maxima <- apply(weights, 1, max)
if (any(rowSums(weights == maxima) != 1)) {
  stop(
    "Each cell must have a unique largest weight to infer its primary cell type.",
    call. = FALSE
  )
}
labels <- setNames(
  colnames(weights)[max.col(weights, ties.method = "first")],
  cells
)

result <- SPLIT::purify(
  counts = counts,
  reference = reference,
  deconvolution_weights = weights,
  primary_cell_type = labels,
  chunk_size = chunk_size,
  DO_output_sce = FALSE
)
purified <- result[["purified_counts"]]
if (
  !setequal(rownames(purified), genes) || !setequal(colnames(purified), cells)
) {
  stop(
    "SPLIT unexpectedly changed the feature or cell identifiers.",
    call. = FALSE
  )
}
purified <- purified[genes, cells, drop = FALSE]
if (any(!is.finite(purified@x)) || any(purified@x < 0)) {
  stop("SPLIT returned invalid purified counts.", call. = FALSE)
}
metadata <- result[["cell_meta"]]
if (!setequal(metadata[["cell_id"]], cells))
  stop("SPLIT returned inconsistent cell metadata.", call. = FALSE)
metadata <- metadata[match(cells, metadata[["cell_id"]]), , drop = FALSE]
outdir <- options[["--outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeMM(purified, file.path(outdir, "matrix.mtx"))
write.table(
  features,
  file.path(outdir, "features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  cells,
  file.path(outdir, "barcodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  metadata,
  file.path(outdir, "cell_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
