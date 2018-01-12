---
title:
- Iris Supplemental Info
author:        
- Brandon Monier
- Adam McDermaid
- Qin Ma
date:       
- 2018-01-12 at 10:36:52
Last Modified: 2018-01-12 at 12:05:17
---

# Accessibility 
IRIS can be freely accessed directly through this [link](http://bmbl.sdstate.edu/VIDGER/) (**change link later**) or through R using the following commands (**not available yet**):

```{r}
if (!require("shiny")) install.packages("shiny")
shiny::runGitHub("vidger-shiny", "btmonier")
```

Typically, the link will provide an easier route to using IRIS. In circumstances where internet connections will be limited (such as during travel), loading IRIS through R while internet is still available will allow users to utilize IRIS without an internet connection later on.

# Input data
IRIS requires two pieces of information for analysis. The first is an expression estimation matrix, also referred to as a count matrix, displaying the gene expression estimates for each sample. The format requires a CSV file with the row names to list the gen IDs and column names to list the sample IDs. The second required input is a condition matrix, wherein the factor levels for each sample are provided. This file requires a CSV format and row names to be the sample IDs matching the sample IDs from the expression estimation matrix and the column names to be the condition factors. 

For the following tutorial, the data used for thhis tutorial are derived from 28 *Vitis vinifera* (i.e. grape) samples with three distinct factors (Rootstock, row, and block).

## Expression matrices
Typically, an expression matrix, also known as count data or a count matrix referes to data where every $i$-th row and $j$-th column refer to how many reads are assigned to gene (ID) $i$ in sample $j$. For example, if we have simplified count data for 4 samples and three genes, the `R` output will look something like this:

```{r}
        sample1 sample2 sample3 sample4
gene001      23       3      45       2
gene002       6       7       7       8
gene003       0      34       3      42
```

**NOTE:** When loading count data into IRIS, make sure that the first column is your gene IDs and that sample names are short, concise, and avoid the use of mathematical operators (`+`, `-`, `/`, `*`, `^`, etc.) and spaces between words. If a space is necessary for legibility, please consider using an underscore (`_`)

## Condition matrices
Condition matrices, also known as metadata,  details the design of your experiment. In this type of matrix, every $i$-th row and $j$-th column refer to factor levels assigned to sample $i$ and factor $j$. For example, if were to look at the samples given in the [count data](#count-mat) section, the metadata `R` output will look something like this:

```{r}
        condition time
sample1   treated   0h
sample2 untreated   0h
sample3   treated  24h
sample4 untreated  24h
```

**NOTE 1:** When loading metadata into IRIS, make sure that the first column is your sample names and that column names and treatment levels are short, concise, and avoid the use of mathematical operators (`+`, `-`, `/`, `*`, `^`, etc.) and spaces between words. If a space is necessary for legibility, please consider using an underscore (`_`) 

**NOTE 2:** Metadata can be expanded to fit the nature of your experiment (i.e. multiple factors can be added). The only thing that must remain consistent between these two matrices, is the sample information. Column names in count data **must** be the same as row names in the metadata.
