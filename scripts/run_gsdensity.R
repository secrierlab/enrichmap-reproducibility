#!/usr/bin/env Rscript
# Run GSDensity (Liang et al. 2023, Nat Commun 14:8255; github.com/qingnanl/gsdensity) on task
# directories exported from Python by revision_utils / gsdensity_bridge.py.
#
# Usage: Rscript run_gsdensity.R <task_dir> [<task_dir> ...]
#
# Each task_dir contains
#   data.mtx      genes x spots matrix of log-normalised expression (Matrix Market)
#   genes.txt     gene names (rows), spots.txt spot names (columns)
#   coords.csv    spot,x,y
#   genesets.json {"name": ["gene", ...], ...}
#   params.json   optional: nn_use (300), restart (0.75), dims (10), spatial (true), spatial_n (10), spatial_ntimes (100), seed (0)
# and receives
#   gsd_scores.csv       per-spot label-propagation (RWR) probabilities, one column per gene set
#   gsd_labels.csv       GSDensity's own positive/negative binarisation (compute.cell.label, multimode)
#   gsd_spatial_kld.csv  pathway-level spatial relevance (compute.spatial.kld: weighted kde2d vs uniform, permutation p)
#   gsd_timing.json      wall-clock seconds per step
#
# The pipeline follows the package defaults and its 10x Visium vignette (compute.mca -> compute.nn.edges(300)
# -> run.rwr.list(restart = 0.75) -> compute.cell.label.df -> compute.spatial.kld(n = 10)). The only change is
# that MCA is run on the supplied log-normalised matrix rather than re-normalising counts inside compute.mca(), so
# that every method in the benchmark sees exactly the same expression values (including noise-perturbed data).

suppressPackageStartupMessages({
  library(gsdensity); library(Seurat); library(CelliD); library(Matrix); library(jsonlite)
})
# dnet with Matrix >= 1.5 may fail on as(<dsCMatrix>, "dgCMatrix"); the README workaround is a dense adjacency
# (igraph.options(sparsematrices = FALSE)). We try sparse first and fall back to dense only on that error.

mca_from_data <- function(X, dims.use = 1:10) {
  # identical to gsdensity::compute.mca() except that NormalizeData() is not re-applied
  obj <- CreateSeuratObject(counts = X)
  obj <- SetAssayData(obj, slot = "data", new.data = X)
  obj <- RunMCA(obj, features = rownames(obj), slot = "data", assay = "RNA")
  g <- intersect(rownames(obj), rownames(obj@reductions$mca@feature.loadings))
  rbind(obj@reductions$mca@feature.loadings[g, dims.use], obj@reductions$mca@cell.embeddings[, dims.use])
}

run_task <- function(td) {
  t <- list(); t0 <- Sys.time()
  p <- list(nn_use = 300, restart = 0.75, dims = 10, spatial = TRUE, spatial_n = 10, spatial_ntimes = 100, seed = 0)
  if (file.exists(file.path(td, "params.json"))) p <- modifyList(p, fromJSON(file.path(td, "params.json")))
  X <- as(readMM(file.path(td, "data.mtx")), "CsparseMatrix")
  genes <- readLines(file.path(td, "genes.txt")); spots <- readLines(file.path(td, "spots.txt"))
  # Seurat replaces '_' by '-' in feature names; apply the same mapping to the gene sets
  genes <- gsub("_", "-", genes)
  rownames(X) <- genes; colnames(X) <- spots
  gs <- lapply(fromJSON(file.path(td, "genesets.json")), function(v) gsub("_", "-", v))
  t$load <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  set.seed(p$seed)
  t1 <- Sys.time(); ce <- mca_from_data(X, dims.use = 1:p$dims); t$mca <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  t1 <- Sys.time(); el <- compute.nn.edges(coembed = ce, nn.use = p$nn_use); t$nn <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  t1 <- Sys.time()
  rwr_all <- function() {
    if (length(gs) > 1) return(run.rwr.list(el = el, gene_set_list = gs, cells = spots, restart = p$restart, parallel = FALSE))
    v <- run.rwr(el = el, gene_set = gs[[1]], cells = spots, restart = p$restart)  # run.rwr.list drops to a vector for 1 set
    setNames(data.frame(v, row.names = spots), names(gs))
  }
  cv <- tryCatch(rwr_all(), error = function(e) {
    message("sparse RWR failed (", conditionMessage(e), "); retrying with dense adjacency")
    igraph::igraph.options(sparsematrices = FALSE); on.exit(igraph::igraph.options(sparsematrices = TRUE))
    rwr_all()
  })
  t$rwr <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  write.csv(cv, file.path(td, "gsd_scores.csv"))

  lab <- sapply(colnames(cv), function(n) tryCatch(compute.cell.label(setNames(cv[[n]], rownames(cv))),
                                                     error = function(e) rep(NA_character_, nrow(cv))))
  rownames(lab) <- rownames(cv)
  write.csv(lab, file.path(td, "gsd_labels.csv"))

  if (isTRUE(p$spatial)) {
    t1 <- Sys.time()
    coords <- read.csv(file.path(td, "coords.csv"), row.names = 1)[rownames(cv), c("x", "y")]
    set.seed(1)
    sk <- t(sapply(colnames(cv), function(n) compute.spatial.kld(spatial.coords = coords, weight_vec = cv[rownames(coords), n],
                                                                 n = p$spatial_n, n.times = p$spatial_ntimes)))
    sk <- as.data.frame(sk); sk$delta.kld <- sk$spatial.kld - sk$rspatial.kld; sk$p.adj <- p.adjust(sk$p.value, "BH")
    write.csv(sk, file.path(td, "gsd_spatial_kld.csv"))
    t$spatial_kld <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  }
  t$total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write(toJSON(t, auto_unbox = TRUE), file.path(td, "gsd_timing.json"))
  cat(sprintf("[gsdensity] %s done in %.1fs (mca %.1fs, rwr %.1fs)\n", td, t$total, t$mca, t$rwr))
}

for (td in commandArgs(trailingOnly = TRUE)) {
  res <- tryCatch(run_task(td), error = function(e) {
    cat(sprintf("[gsdensity] %s FAILED: %s\n", td, conditionMessage(e)))
    writeLines(conditionMessage(e), file.path(td, "gsd_error.txt"))
  })
}
