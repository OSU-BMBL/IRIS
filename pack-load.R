#---------------------------------------------------------------------
# Title:         IRIS - Package Loader
# Author:        Brandon Monier
# Created:       2018-03-14 at 17:02:31
# Last Modified: 2018-11-07 at 10:35:54
#---------------------------------------------------------------------

packages <- c(
	"crosstalk", "dplyr", "DT", "gtools", "plotly", "shiny", "plyr",
	"shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
	"Rcpp", "Hmisc", "ggplot2", "geneplotter", "locfit", "GGally", 
	"pheatmap",	"reshape2", "backports", "digest", "fields", "psych",
	"DESeq2", "edgeR", "limma", "QUBIC", "stringr", "tools", "openxlsx",
	"Rtsne", "WGCNA", "flashClust", "parallel", "MCL", "kmed"
)

lapply(packages, require, character.only = TRUE)
