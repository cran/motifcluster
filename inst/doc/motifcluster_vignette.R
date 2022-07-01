## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#> "
)
set.seed(200401293)
old_options <- options()
options(digits = 2)
library(mclust)
# nolint start: line_length_linter

## ---- eval = FALSE------------------------------------------------------------
#  library(devtools)
#  library(mclust)

## ---- eval = FALSE------------------------------------------------------------
#  install_github("wgunderwood/motifcluster/R")

## ---- results = "hide"--------------------------------------------------------
library(motifcluster)

## -----------------------------------------------------------------------------
G1 <- matrix(c(
  0, 2, 0, 0,
  0, 0, 2, 3,
  0, 4, 0, 0,
  4, 0, 5, 0
), byrow = TRUE, nrow = 4)

## -----------------------------------------------------------------------------
get_motif_names()

## -----------------------------------------------------------------------------
build_motif_adjacency_matrix(G1, motif_name = "M1")

## -----------------------------------------------------------------------------
build_motif_adjacency_matrix(G1, motif_name = "M1", motif_type = "func")

## -----------------------------------------------------------------------------
build_motif_adjacency_matrix(G1, motif_name = "M1", motif_type = "func",
  mam_weight_type = "mean")

## -----------------------------------------------------------------------------
build_motif_adjacency_matrix(G1, motif_name = "M1", motif_type = "func",
  mam_weight_type = "product")

## -----------------------------------------------------------------------------
mam_sparse <- build_motif_adjacency_matrix(G1, motif_name = "M1", mam_method = "sparse")
mam_dense  <- build_motif_adjacency_matrix(G1, motif_name = "M1", mam_method = "dense")
all(mam_sparse == mam_dense)

## -----------------------------------------------------------------------------
block_sizes <- c(5, 3)
connection_matrix <- matrix(c(
  0.9, 0.2,
  0.2, 0.9
), nrow = 2, byrow = TRUE)

## -----------------------------------------------------------------------------
sample_dsbm(block_sizes, connection_matrix)

## -----------------------------------------------------------------------------
weight_matrix <- matrix(c(
  5, 2,
  2, 5
), nrow = 2, byrow = TRUE)
sample_dsbm(block_sizes, connection_matrix, weight_matrix,
  sample_weight_type = "constant")

## -----------------------------------------------------------------------------
sample_dsbm(block_sizes, connection_matrix, weight_matrix,
  sample_weight_type = "poisson")

## -----------------------------------------------------------------------------
source_block_sizes <- c(2)
destination_block_sizes <- c(3, 2)
bipartite_connection_matrix <- matrix(c(0.9, 0.2), nrow = 1)
sample_bsbm(source_block_sizes, destination_block_sizes, bipartite_connection_matrix)

## -----------------------------------------------------------------------------
bipartite_weight_matrix <- matrix(c(7, 2), nrow = 1)
sample_bsbm(source_block_sizes, destination_block_sizes, bipartite_connection_matrix,
  bipartite_weight_matrix, sample_weight_type = "poisson")

## -----------------------------------------------------------------------------
G2 <- matrix(c(
  0, 2, 0, 0,
  2, 0, 4, 3,
  0, 4, 0, 5,
  0, 3, 5, 0
), byrow = TRUE, nrow = 4)

## -----------------------------------------------------------------------------
build_laplacian(G2, type_lap = "comb")

## -----------------------------------------------------------------------------
build_laplacian(G2, type_lap = "rw")

## -----------------------------------------------------------------------------
spectrum <- run_laplace_embedding(G2, num_eigs = 2, type_lap = "rw")
spectrum$vals
spectrum$vects

## -----------------------------------------------------------------------------
G3 <- matrix(c(
  0, 0, 0, 0,
  0, 0, 2, 3,
  0, 4, 0, 0,
  4, 0, 5, 0
), byrow = TRUE, nrow = 4)

## -----------------------------------------------------------------------------
spectrum <- run_motif_embedding(G3, motif_name = "M1", motif_type = "func",
  mam_weight_type = "unweighted", mam_method = "sparse", num_eigs = 2,
  restrict = TRUE, type_lap = "rw")
spectrum$vals
spectrum$vects

## -----------------------------------------------------------------------------
block_sizes <- rep(10, 3)

## -----------------------------------------------------------------------------
connection_matrix <- matrix(c(
  0.8, 0.2, 0.2,
  0.2, 0.8, 0.2,
  0.2, 0.2, 0.8
), nrow = 3)

## -----------------------------------------------------------------------------
weight_matrix <- matrix(c(
  20, 10, 10,
  10, 20, 10,
  10, 10, 20
), nrow = 3)

G4 <- sample_dsbm(block_sizes, connection_matrix, weight_matrix,
  sample_weight_type = "poisson")

## -----------------------------------------------------------------------------
motif_cluster <- run_motif_clustering(G4, motif_name = "M1", motif_type = "func",
  mam_weight_type = "mean", mam_method = "sparse", type_lap = "rw",
  num_eigs = 4, num_clusts = 3
)

## -----------------------------------------------------------------------------
truth <- c(rep(1, 10), rep(2, 10), rep(3, 10))
mclust::adjustedRandIndex(motif_cluster$clusts, truth)

## ---- include = FALSE---------------------------------------------------------
options(old_options)
#nolint end

