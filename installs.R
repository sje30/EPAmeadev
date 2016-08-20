## install steps for this code.

install.packages(c("doBy", "ggplot2", "gridExtra", "xtable"))
install.packages(c("pracma", "randomForest", "e1071"))


## These two packages take a while to install, so do it
## only if they are missing.

if (!require(rhdf5)) {
  source("http://bioconductor.org/biocLite.R"); biocLite("rhdf5")
}

if (!require(sjemea)) {
  devtools::install_github("sje30/sjemea",build_vignettes=TRUE)
}


