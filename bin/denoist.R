#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(DenoIST)
    library(SpatialExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, required = TRUE, default = NULL) {
    i <- match(name, args)
    if (is.na(i) || i == length(args)) {
        if (required) stop("Missing required argument: ", name)
        return(default)
    }
    args[[i + 1]]
}

mat_path <- get_arg("--mat")
tx <- readRDS(get_arg("--tx"))
coords_path <- get_arg("--coords", required = FALSE)

mat <- readRDS(mat_path)
if (inherits(mat, "SpatialExperiment")) {
    # DenoIST reads counts and spatial coordinates from the SpatialExperiment.
    coords <- NULL
} else {
    if (is.null(coords_path)) stop("--coords is required when --mat is a count matrix.")
    coords <- readRDS(coords_path)
}

result <- denoist(
    mat = mat,
    tx = tx,
    coords = coords,
    tx_x = get_arg("--tx-x", default = "x"),
    tx_y = get_arg("--tx-y", default = "y"),
    distance = 50,
    nbins = 200,
    cl = 1,
    out_dir = "denoist_results",
    neighbour_mode = "fast"
)

saveRDS(result, get_arg("--output"))
