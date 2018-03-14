---
 title:         IRIS README
 author:        Brandon Monier
 created:       2018-01-12 at 09:56:43
 last modified: 2018-03-14 at 17:34:12
---

# IRIS
<http://bmbl.sdstate.edu/IRIS/>

## About
IRIS-DGE (**I**nteractive **R**NA-seq analysis and **I**nterpretation using
**S**hiny-**D**ifferential **G**ene **E**xpression), is a web-based tool for
the analysis of RNA-seq count data. This tool’s purpose is to provide users 
with a comprehensive and user-friendly method for performing differential gene 
expression (DGE) analysis regardless of their computational experience. 
IRIS-DGE also has integrated experimental design options to cater to users 
with non-traditional DGE requirements, such as interaction terms or paired 
data. This tool is designed in a way for usable results to be generated in 
around one minute or for users to invest more time into detailed 
investigations of their data. IRIS is a **user-friendly** and **interactive** 
Shiny app for gene expression analysis. This app takes advantage of several 
popular DGE tools (*DESeq2*, *edgeR*, and *limma*) available through 
Bioconductor in conjunction with the Plotly and DataTable API libraries for R.

## Get local application
To get a local version of IRIS-DGE, simply copy and paste this code into an R
terminal:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("iris", "btmonier")
```