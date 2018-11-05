#---------------------------------------------------------------------
# Title:         IRIS - Package Loader
# Author:        Brandon Monier
# Created:       2018-03-14 17:02:31 CDT
# Last Modified: 2018-05-22 15:26:55 CDT 
#---------------------------------------------------------------------

packages <- c(
	"crosstalk", "dplyr", "DT", "gtools", "plotly", "shiny", "plyr",
	"shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
	"Rcpp", "Hmisc", "ggplot2", "geneplotter", "locfit", "GGally", 
	"pheatmap",	"reshape2", "backports", "digest", "fields", "psych",
	"DESeq2", "edgeR", "limma", "QUBIC", "stringr", "tools", "openxlsx",
	"Rtsne", "WGCNA", "flashClust"
)

lapply(packages, require, character.only = TRUE)