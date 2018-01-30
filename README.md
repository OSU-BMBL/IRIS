---
 Title:         IRIS Documentation
 Author:        Brandon Monier
 Created:       2018-01-12 at 09:56:43
 Last Modified: 2018-01-30 at 15:20:04
---

# IRIS

## Link
[bmbl.sdstate.edu/IRIS/](http://bmbl.sdstate.edu/IRIS/)

## About
IRIS (**I**nteractive **R**NA-seq analysis and **I**nterpretation using
**S**hiny), is a web-based tool for the analysis of RNA-seq count data. This
tool’s purpose is to provide users with a comprehensive and user-friendly
method for performing differential gene expression (DGE) analysis regardless of
their computational experience. IRIS also has integrated experimental design
options to cater to users with non-traditional DGE requirements, such as
interaction terms or paired data. This tool is designed in a way for usable
results to be generated in around one minute or for users to invest more time
into detailed investigations of their data. IRIS is a **user-friendly** and
**interactive** Shiny app for gene expression analysis. This app takes
advantage of several popular DGE tools (*DESeq2*, *edgeR*, and *limma*)
available through Bioconductor in conjunction with the Plotly and DataTable API
libraries for R. 

## Local installation
To run the the application locally, you can install the `shiny` package in
**R**, and use the function `runGithub()`:

``` r
# Placeholder only!
if (!require("shiny")) install.packages("shiny")
shiny::runGitHub("iris", "btmonier")
```

## To-do
- [ ] Finish supplemental information (`/vignette/iris-supplemental.pdf`)
- [ ] Transfer source code to new directory
- [ ] Finish paper
