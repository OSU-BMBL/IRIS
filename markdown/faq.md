## **F**requently **A**sked **Q**uestions

#### Table of Contents
* [What does VIDGER stand for?](#vidname)
* [What web browser(s) can I use?](#web-browser)
* [What does DEG/DGE mean?](#deg)
* [What methods does VIDGER use for differential gene expression?](#dgemethods)
* [What is count data?](#count-mat)
* [What is metadata?](#metadata)
* [What do you mean by "filter cutoff"?](#filtercutoff)
* [On the 'Submit and QC' page, what are the transformation methods for counts?](#subqctran)
* [What is a correlation matrix?](#cormatrix)
* [What is a sample distance matrix?](#sampdistmatrix)
* [What is PCA?](#vis-pca)
* [What is MDS?](#vis-mds)
* [What are heatmaps?](#heatmaps)
* [What is Biclustering?](#biclustering)
* [What experimental designs does VIDGER currently allow?](#expdesign)
* [What is a Adj. $p$-value cutoff?](#adjpval)
* [What is a Min. fold change value?](#minlfc)
* [What is an MA plot?](#maplot)
* [What is a volcano plot?](#volplot)
* [Where can I download a local version of the app?](#localapp)
* [Where can I download the ViDGER package?](#vidgerpackage)
* [References](#refs)

***
#### What does VIDGER stand for? <a id="vidname"></a>
VIDGER stands for **V**isualization and **I**nterpretation of **D**ifferential **G**ene **E**xpression using **R**. This web application is based off of our prior package [ViDGER](https://github.com/btmonier/vidger) which stands for **Vi**sualization of **D**ifferential **G**ene **E**xpression using **R**

<br>

#### What web browser(s) can I use? <a id="web-browser"></a>
VIDGER has been properly tested on both [Firefox](https://www.mozilla.org/en-US/firefox/) and [Chrome](https://www.google.com/chrome/browser/desktop/index.html), so we recommend using either of these browsers.

<br>

#### What does DEG/DGE mean? <a id="deg"></a>
DEG and DGE simply mean **d**ifferential **e**xpression of **g**enes, and **d**ifferential **g**ene **e**xpression, respectively.

<br>

#### What methods does VIDGER use for differential gene expression? <a id="dgemethods"></a>
VIDGER currently uses [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), and [limma](http://bioconductor.org/packages/release/bioc/html/limma.html) (Love et al. 2014; McCarthy et al. 2012; Ritchie et al. 2015; Robinson et al. 2010)

<br>

#### What is count data? <a id="count-mat"></a>
Count data, also referred to as an expression or count matrix refers to data where every $i$-th row and $j$-th column refer to how many reads are assigned to gene (ID) $i$ in sample $j$. For example, if we have simplified count data for 4 samples and three genes, the `R` output will look something like this:

```
        sample1 sample2 sample3 sample4
gene001      23       3      45       2
gene002       6       7       7       8
gene003       0      34       3      42
```

**NOTE:** When loading count data into VIDGER, make sure that the first column is your gene IDs and that sample names are short, concise, and avoid the use of mathematical operators (`+`, `-`, `/`, `*`, `^`, etc.) and spaces between words. If a space is necessary for legibility, please consider using an underscore (`_`) 

<br>

#### What is metadata? <a id="metadata"></a>
Metadata, also known as a condition matrix, details the design of your experiment. In this type of matrix, every $i$-th row and $j$-th column refer to factor levels assigned to sample $i$ and factor $j$. For example, if were to look at the samples given in the [count data](#count-mat) section, the metadata `R` output will look something like this:

```
        condition time
sample1   treated   0h
sample2 untreated   0h
sample3   treated  24h
sample4 untreated  24h
```

**NOTE 1:** When loading metadata into VIDGER, make sure that the first column is your sample names and that column names and treatment levels are short, concise, and avoid the use of mathematical operators (`+`, `-`, `/`, `*`, `^`, etc.) and spaces between words. If a space is necessary for legibility, please consider using an underscore (`_`) 

**NOTE 2:** Metadata can be expanded to fit the nature of your experiment (i.e. multiple factors can be added). The only thing that must remain consistent between these two matrices, is the sample information. Column names in count data **must** be the same as row names in the metadata.

<br>

#### What do you mean by "filter cutoff"? <a id="filtercutoff"></a>
A "filter cutoff" is a numerical parameter ($n$) will be implemented in your count data which will be used to remove any genes or ID rows that are $< n$. The default parameter is set to $10$, but this can be changed to any integer $>$ $0$. If you wish to avoid this step, simply set the filter cutoff to $0$.

<br>

#### On the 'Submit and QC' page, what are the transformation methods for counts? <a id="subqctran"></a>
VIDGER can transform count using the transformation functions from the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package (Love et al. 2014). Currently there are four options:

* **Normal log: log2(n + pseudocount)**  
  The most common transformation technique is the logarithm, defined as "Normal log" in the parameters. This is implemented using the following formula:  
  $$y = log_{2}(n + p)$$
  Where $n$ is the read counts and $p$ is a "pseudocount" due to some gene ID read counts being $0$.

* **Regularized log: rlog(n)**  
	Genes that posses high levels of raw counts will give similar results to  normalization techniques, however, for genes that have lower levels of counts, the overall values are reduced to fit more to the average values produced by genes across entire samples in an experimental setup. The overall process is defined by the following formula:  
	$$log_{2}(q_{ij}) = \beta_{i0} + \beta_{ij}$$
	Where $q_{ij}$ is relative to the expected value of counts for gene $i$ and sample $j$, $\beta_{i0}$ is the intercept, $\beta_{ij}$ is the shrunken effect for each sample based on a dispersion-mean trend used on the count data.

	**NOTE:** the `rlog()` function is rather time consuming, so be prepared to wait a while depending on the size of your experiment.

* **Variance stabilizing transform: vst(n)**     
	Variance stabilizing transformation (VST) is calculated from the fitted dispersion-mean relations and then transforms the normalized count data which produces homoskedastic data. The overall goal in using VST is to produce constant variance that is relative to the mean (Tibshirani 1988).

* **No transformation**  
	If you choose this option, no transformation will applied to your count data.

These transformation options are not directly applied to the data used in the "DGE Analysis" tab. These are used for visualization and biclustering purposes in the "Preliminary Analysis" tab and count summary visualizations.

<br>

#### What is a correlation matrix? <a id="cormatrix"></a>
A correlation matrix is used to generate the correlation coefficients between all relationships amongst the samples in your data. To visualize this data, a heatmap is used. The intense the yellow color is in the matrix refers to a higher positive correlation between two samples while the more blue refers to lower coefficients.

**NOTE:** This matrix utilizes the Plotly API library for R, so it is interactive. Click on the cells of the matrix to generate an interactive scatterplot for the two relative samples.

<br>

#### What is a sample distance matrix? <a id="sampdistmatrix"></a>
Sample correlation matrices are used to determine Euclidean distances between the gene expression vectors for each sample pair. To visualize these relationships, a heatmap is used. The more blue a cell is indicates that the sample pairs overall gene expression is the same. The more intense the red color becomes indicates that gene expression is widely different. Additionally, sample distances are clustered to show additional groupings amongst samples.

<br>

#### What is PCA? <a id="vis-pca"></a>
PCA stands for **P**rincipal **C**omponent **A**nalysis. This statistical technique is used to reduce complex datasets of multiple variables of gene expression to fewer dimensions (plotted on $x$, $y$, or $z$ axes). Generally, for RNA-seq data, it is used to see clustering or grouping amongst sample and their respective treatments. For example, if you have a "treated" and "untreated" conditional levels, you would generally expect samples from each of the levels to group closer together.

<br>

#### What is MDS? <a id="vis-mds"></a>
Similar to PCA, **m**ultiple-**d**imensional **s**caling (MDS) is another grouping technique for samples. Unlike PCA, MDS produces components through non-linear methods for gene expression.

<br>

#### What are heatmaps? <a id="heatmaps"></a>
Heatmaps, specifically under the "Heatmap" tab under "Preliminary Analysis" are used to visualize the gene IDs that have the highest overall mean (i.e. row means in the count data) amongst the samples in a RNA-seq experiment. The more intense the yellow color indicates higher transformed counts (see the "[On the 'Submit and QC' page, what are the transformation methods for counts?](#subqctran)" question for further information). Since count data matrices can become incredibly large in size, showing the entire data set would too "noisy". As a result, we have set an ID cutoff parameter to show the $20$ most expressed IDs for display. This parameter can be changed to any integer $> 0$. 

Since this interactive, you can hover over cells for additional information and click on cells to generate an interactive dot plot which will display total read counts for that specific gene ID by the selected factor.

<br>

#### What is biclustering? <a id="biclustering"></a>
Biclustering is a clustering technique which allows for clustering of both gene IDs and samples in a count data matrix. Biclustering in VIDGER performs a biclustering analysis using one of the selected biclustering tools with a maximum bicluster size of the indicated cutoff value. 

<br>

#### What experimental designs does VIDGER currently allow? <a id="expdesign"></a>
One of the motivations when designing VIDGER was to allow for a variety of commonly used experimental designs employed in RNA-seq experiments. Currently, our application allows for:

* Comparisons between two factor levels (groups) amongst one treatment factor; 
* Comparisons between level combinations between two treatment factors; 
* Classical interaction designs between two treatment factors; 
* Accountability for additive models (blocking and paired effects)

<br>

#### What is a Adj. $p$-value cutoff? <a id="adjpval"></a>
The adjusted (adj.) $p$-value cutoff is a numerical parameter to filter gene IDs that are not significantly expressed between two conditions. Any adjusted $p$-value that is $< n$, where $n$ is the cutoff parameter, will be removed from the filtered data and visualizations in conjunction with the [minimum $log_{2}$ fold change](#minlfc).

<br>

#### What is a Min. fold change value? <a id="minlfc"></a>
The minimum (Min.) fold change is a numerical parameter to filter gene IDs that have a $log_{2}$ fold change that is $< n$, where $n$ is the cutoff parameter. Any gene ID that is below the parameter will be removed from the filtered data and visualizations in conjunction with the [adjusted $p$-value](#adjpval).

<br>

#### What is an MA plot? <a id="maplot"></a>
An MA plot is a method to visualize gene expression data between two conditions. A scatter plot is produced by plotting the $log_{2}$ fold change (**M**) compared to the transformed base mean ($log_{10}$) (**A**) on the $x$ axis.

<br>

#### What is a volcano plot <a id="volplot"></a>
A volcano plot is another method to visualize gene expression data between two conditions. In this plot, the transformed $p$-values ($-log_{10}$) is compared to the $log_{2}$ fold change. Volcano plots can be used for rapid identification of statistically significant genes in terms of there $p$-value and $log_{2}$ fold change. Significant values will generally be located in the upper right- and left-handed corners of the graph. 

<br>

#### Where can I download a local version of the app? <a id="localapp"></a>

**GitHub**

You can download the latest (*and experimental*) version of the web application using this script in an up-to-date version of `R`:

```{r}
if (!require("shiny")) install.packages("shiny")
shiny::runGitHub("vidger-shiny", "btmonier")
```

<br>

#### Where can I download the ViDGER package? <a id="vidgerpackage"></a>

**GitHub**

You can download the latest (*and experimental*) GitHub repo this web application was based on using this script in an up-to-date version of `R`:

```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("btmonier/vidger")
```

<br>

#### References <a id="refs"></a>
Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology 15 (12): 550.

McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), pp. 4288-4297.

Tibshirani, Robert. 1988. “Estimating Transformations for Regression via Additivity and Variance Stabilization.” Journal of the American Statistical Association 83: 394–405.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), pp. e47.

Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), pp. 139-140.

<br>
<br>
<br>
