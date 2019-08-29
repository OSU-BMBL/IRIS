#---------------------------------------------------------------------
# Title:         IRIS - Package Loader
# Author:        Brandon Monier
# Created:       2018-03-14 at 17:02:31
# Last Modified: 2019-08-29 at 10:01:05
#---------------------------------------------------------------------

packages <- c(
    "crosstalk", "plyr", "DT", "gtools", "plotly", "shiny", "dplyr",
    "shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
    "Rcpp", "Hmisc", "ggplot2", "geneplotter", "locfit", "GGally",
    "pheatmap", "reshape2", "backports", "digest", "fields", "psych",
    "DESeq2", "edgeR", "limma", "QUBIC", "stringr", "tools", "openxlsx",
    "Rtsne", "WGCNA", "flashClust", "parallel", "MCL", "kmed", "ape",
    "BRIC", "GO.db"
)

lapply(packages, require, character.only = TRUE)
