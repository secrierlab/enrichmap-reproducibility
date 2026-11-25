#!/bin/bash
# Build an isolated R environment for GSDensity (Liang et al. 2023, Nat Commun; github.com/qingnanl/gsdensity).
# GSDensity is an R package; dnet and anticlust (<=0.6.1) are archived on CRAN and are installed from the archive.
set -e
eval "$(conda shell.bash hook)"
mamba create -y -n gsdensity_env -c conda-forge -c bioconda \
  r-base=4.3 r-seurat=4.4 r-seuratobject=4.1 bioconductor-cellid bioconductor-rgraphviz bioconductor-suprahex bioconductor-graph \
  r-igraph r-infotheo r-multimode r-philentropy r-rann r-future.apply r-matrix r-dplyr r-reshape2 r-ggplot2 \
  r-remotes r-data.table r-ape r-hexbin r-jsonlite r-irkernel r-mass
conda activate gsdensity_env
Rscript -e "options(repos=c(CRAN=\"https://cloud.r-project.org\")); Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=\"true\"); \
  remotes::install_version(\"dnet\", version=\"1.1.7\", upgrade=\"never\", dependencies=FALSE); \
  remotes::install_version(\"anticlust\", version=\"0.6.1\", upgrade=\"never\", dependencies=TRUE); \
  remotes::install_github(\"qingnanl/gsdensity\", upgrade=\"never\", dependencies=FALSE); \
  library(gsdensity); cat(\"GSDensity installed:\", as.character(packageVersion(\"gsdensity\")), \"\\n\")"
