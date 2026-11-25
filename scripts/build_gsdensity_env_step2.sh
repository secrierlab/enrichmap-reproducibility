#!/bin/bash
# Step 2: archived CRAN deps + GSDensity, installed strictly into the env library (see activate.d/zz_isolate_r.sh).
# anticlust is installed with hard dependencies only (its Suggests, e.g. testthat/Rglpk, are not needed).
set -e
eval "$(conda shell.bash hook)"
conda activate gsdensity_env
Rscript -e "L <- .libPaths()[1]; stopifnot(grepl(\"gsdensity_env\", L)); options(repos=c(CRAN=\"https://cloud.r-project.org\")); Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=\"true\"); \
  hard <- c(\"Depends\", \"Imports\", \"LinkingTo\"); \
  remotes::install_version(\"dnet\", version=\"1.1.7\", upgrade=\"never\", dependencies=hard, lib=L); \
  remotes::install_version(\"anticlust\", version=\"0.6.1\", upgrade=\"never\", dependencies=hard, lib=L); \
  remotes::install_github(\"qingnanl/gsdensity\", upgrade=\"never\", dependencies=FALSE, lib=L); \
  library(gsdensity); cat(\"GSDensity installed:\", as.character(packageVersion(\"gsdensity\")), \"\\n\"); print(sessionInfo())"
