
IRIS <img src="www/logo.svg" align="right" height="120"/>
=========================================================

<https://bmbls.bmi.osumc.edu/IRIS/>

Alternate server: <http://bmbl.sdstate.edu/IRIS/>

Citation
--------

### APA

Monier, B., McDermaid, A., Wang, C., Zhao, J., Miller, A., Fennell, A., & Ma, Q. (2019). [IRIS-EDA: an integrated RNA-Seq interpretation system for gene expression data analysis.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006792) *PLoS computational biology*, 15(2), e1006792.

### BibTeX

    @article{monier2019iris,
      title={IRIS-EDA: an integrated RNA-Seq interpretation system for gene expression data analysis},
      author={Monier, Brandon and McDermaid, Adam and Wang, Cankun and Zhao, Jing and Miller, Allison and Fennell, Anne and Ma, Qin},
      journal={PLoS computational biology},
      volume={15},
      number={2},
      pages={e1006792},
      year={2019},
      publisher={Public Library of Science}
    }

Overview
--------

IRIS-EDA (**I**nteractive **R**NA-seq analysis and **I**nterpretation using **S**hiny-**E**xpression **D**ata **A**nalysis), is a web-based tool for the analysis of RNA-seq count data. This tool's purpose is to provide users with a comprehensive and user-friendly method for performing differential gene expression (DGE) analysis regardless of their computational experience. IRIS-EDA also has integrated experimental design options to cater to users with non-traditional DGE requirements, such as interaction terms or paired data. This tool is designed in a way for usable results to be generated in around one minute or for users to invest more time into detailed investigations of their data. IRIS is a **user-friendly** and **interactive** Shiny app for gene expression analysis. This app takes advantage of several popular DGE tools (*DESeq2*, *edgeR*, and *limma*) available through Bioconductor in conjunction with the Plotly and DataTable API libraries for R.

Installation
------------

To get a **local version** of IRIS-EDA, simply copy and paste the following code chunks into an R terminal:

### Install CRAN packages

IRIS-EDA requires several packages to operate. Run this code to get the necessary packages from the CRAN repository:

``` r
# CRAN
packages <- c(
    "crosstalk", "dplyr", "DT", "gtools", "plotly", "shiny", "plyr",
    "shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
    "Rcpp", "Hmisc", "ggplot2", "locfit", "GGally", "pheatmap",
    "reshape2", "backports", "digest", "fields", "psych", "stringr",
    "tools", "openxlsx", "Rtsne", "WGCNA", "flashClust", "parallel",
    "MCL", "kmed", "ape"
)
np <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(np)) install.packages(np)
```

### Install Bioconductor packages

You will also need several Bioconductor packages. Run this code to get the necessary packages from the Bioconductor repository:

``` r
# Bioconductor
bioc_packages <- c(
    "DESeq2", "edgeR", "limma", "QUBIC", "geneplotter", "GO.db", "impute",
    "preprocessCore", "AnnotationDbi"
)
np <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(np)
```

### Install Developmental packages

To run BRIC analysis, you also need to download the source code for this clustering algorithm. Run this code to get the GitHub package:

``` r
# GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("OSU-BMBL/BRIC", force = T)
```

### Run the Shiny application

Once you have installed all of the necessary packages, you can run this code to launch the Shiny application:

``` r
shiny::runGitHub("iris", "OSU-BMBL")
```

------------------------------------------------------------------------

*Last updated:* 2019-08-29
