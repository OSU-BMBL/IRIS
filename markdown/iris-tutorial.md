## IRIS-EDA Tutorial

### Table of Contents
1. [Accessibility](#accessibility)
2. [A note about input data](#input-data)
	* [Expression matrices](#id-emat)
	* [Condition matrices](#id-conmat)
3. [Expedited Analysis](#expedited-analysis)
	* [Overview](#ea-over)
	* [Submit and QC](#ea-subqc)
	* [Data export](#ea-dataexp)
4. [In-depth Analysis](#indepth-analysis)
	* [Overview](#ia-over)
	* [Submit and QC](#ia-subqc)
		* [File Summary](#ia-file)
		* [Count Summary](#ia-count)
	* [Discovery-Driven Analyses](#ia-prelim-analysis)
	* [DGE Analysis](#ia-dge-analysis)
		* [An overview of experimental designs](#ia-over-exp)
		* [Setting up an experiment](#ia-set-exp)
5. [GEO usage](#geo-main)
  * [Overview](#geo-over)
  * [Series](#geo-series)
  * [Samples](#geo-samples)
  * [Protocols](#geo-protocols)
  * [Data processing pipeline](#geo-pipe)
  * [Processed data files](#geo-proc-data)
  * [Raw files](#geo-raw-data)
  * [Paired-end experiments](#geo-pair)
  * [Downloads](#geo-downloads)
6. [Single-cell RNA-seq](#scrna-main)
  * [Overview](#scrna-over)
  * [A method to obtain this data](#scrna-method-01)
  * [TPM filtration](#scrna-tpm-filter)

- - -

### Accessibility <a id="accessibility"></a>
IRIS-EDA can be freely accessed directly through this
[link](http://bmbl.sdstate.edu/IRIS/) or through R using the following commands:

```{r, eval=FALSE}
if (!require("shiny")) install.packages("shiny")
shiny::runGitHub("iris", "btmonier")
```

Typically, the link will provide an easier route to using IRIS-EDA. In
circumstances where internet connections will be limited (such as during
travel), loading IRIS-EDA through R while internet is still available will allow
users to utilize IRIS-EDA without an internet connection later on.

<br/>

### A note about input data <a id="input-data"></a>
IRIS-EDA requires two pieces of information for analysis. The first is an
expression estimation matrix, also referred to as a count matrix, displaying
the gene expression estimates for each sample. The format requires a CSV file
with the row names to list the gen IDs and column names to list the sample IDs.
The second required input is a condition matrix, wherein the factor levels for
each sample are provided. This file requires a CSV format and row names to be
the sample IDs matching the sample IDs from the expression estimation matrix
and the column names to be the condition factors.

For the following tutorial, the data used for this tutorial are derived from
28 *Vitis vinifera* (grape) samples with three distinct factors
(Rootstock, row, and block). This data can viewed as the "big" example data set
found under the `Submit and QC` tab under `1. Submission Parameters`.

<br/>

#### Expression matrices <a id="id-emat"></a>
Typically, an expression matrix, also known as count data or a count matrix
refers to data where every $i$-th row and $j$-th column refer to how many
reads are assigned to gene (ID) $i$ in sample $j$. For example, if we have
simplified count data for 4 samples and three genes, the `R` output will look
something like this:

```
        sample1 sample2 sample3 sample4
gene001      23       3      45       2
gene002       6       7       7       8
gene003       0      34       3      42
...         ...     ...     ...     ...
```

**Note:** When loading count data into IRIS-EDA, make sure that the first column is
your gene IDs and that sample names are short, concise, and avoid the use of
mathematical operators (`+`, `-`, `/`, `*`, `^`, etc.) and spaces between
words. If a space is necessary for legibility, please consider using an
underscore (`_`)

<br/>

#### Condition matrices <a id="id-conmat"></a>
Condition matrices, also known as metadata,  details the design of your
experiment. In this type of matrix, every $i$-th row and $j$-th column refer to
factor levels assigned to sample $i$ and factor $j$. For example, if were to
look at the samples given in the [count data](#count-mat) section, the metadata
`R` output will look something like this:

```
        condition time
sample1   treated   0h
sample2 untreated   0h
sample3   treated  24h
sample4 untreated  24h
```

**Note^1^:** When loading metadata into IRIS-EDA, make sure that the first column is
your sample names and that column names and treatment levels are short,
concise, and avoid the use of mathematical operators (`+`, `-`, `/`, `*`, `^`,
etc.) and spaces between words. If a space is necessary for legibility, please
consider using an underscore (`_`)

**Note^2^:** Metadata can be expanded to fit the nature of your experiment (i.e.
multiple factors can be added). The only thing that must remain consistent
between these two matrices, is the sample information. Column names in count
data **must** be the same as row names in the metadata.

<br/>

<br/>

### Expedited analysis <a id="expedited-analysis"></a>

#### Overview <a id="ea-over"></a>
Expedited analysis is for users who want a quick and efficient method of
producing DGE results using the default parameters and tools in IRIS-EDA.



#### Submit and QC <a id="ea-subqc"></a>
1. Click on the `Submit and QC` tab:

<img src="../vignettes/img/subqc-01.png"/>

2. Load user data by selecting `Load my own data` under
   `1. Submission Parameters`. User data requires one count matrix and one
   condition matrix. **Note:** you must supply both forms of data:

<img src="../vignettes/img/subqc-02.png"/>

3. Click `Submit` to load the data under `3. Launch Overview`:

<img src="../vignettes/img/subqc-03.png"/>

4. After submitting the data, proceed to the `DGE Analysis` tab at the top.
   **Note:** if a specific experimental design other than the basic two-group
   comparison is required, select it accordingly:

<img src="../vignettes/img/subqc-04.png"/>

5. Select the factor of interest and check the boxes for the comparisons of
   interest under `3. Experimental Parameters`. **Note^1^:** you must choose at
   least one comparison for whatever experimental design you choose under
   `1. Experimental Setup`! **Note^2^:** If you do not want to use `DESeq2`,
   feel free to change this method to either `edgeR` or `limma-voom` under
   `2. DGE Paramters`. Also, feel free to change adjusted $p$-values and
   minimum fold change if these defaults do not fit your analysis:

<img src="../vignettes/img/subqc-05.png"/>

5. Submit the parameters to perform DGE analysis under `4. Launch Analysis`:

<img src="../vignettes/img/subqc-03.png"/>

6. Select the newly formed `Plots` subtab on the main page to view the DGE
   results table towards the bottom:

<br/>

<img src="../vignettes/img/subqc-06.png"/>

<br/>

#### Data export <a id="ea-dataexp"></a>
1. To export results data from the `DGE Analysis` tab, click on the
   `Download All Data` button at the bottom of the page to download the
   results file:

<img src="../vignettes/img/subqc-07.png"/>

<br/>

<br/>

### In-depth analysis <a id="indepth-analysis"></a>

#### Overview <a id="ia-over"></a>
By choosing the in-depth analysis route, you will be given more information
about your data including:

* Count data distribution;
* Total number of reads/sample;
* Sample correlation analysis;
* Biclustering;
* Principal components and multidimensional scaling;
* Identification of most variable IDs in data.

<br/>

#### Submit and QC <a id="ia-subqc"></a>
1. Click on the `Submit and QC` tab near the top-left corner of the application:

<img src="../vignettes/img/subqc-01.png"/>

2. Under `1. Submission Paramters`, select either
   `Start with a small example data set`, `Start with a big example data set`,
   or `Load my own data` to upload user data. **Note:** User data requires one
   count matrix and one condition matrix:

<img src="../vignettes/img/subqc-02.png"/>

2. Under `2. Data Processing`, select a filter cutoff to simplify and
   expedite computations by eliminating any rows that have below specified
   expression totals. The default argument for our application is `10`:

<img src="../vignettes/img/subqc-08.png"/>

3. Right below the filter cutoff parameter, select a transformation method for
   the count data for displaying the expression estimate. **Note:** For more
   information about any of these topics, take a look at the *FAQ* section at
   the bottom of this document:

<img src="../vignettes/img/subqc-09.png"/>

4. Click `Submit` to load the data under `3. Launch Overview`:

<img src="../vignettes/img/subqc-03.png"/>

5. After you click `Submit`, the main page will now be populated with several
   pieces of information in two sub-tabs, `File Summary` and `Count Summary`.

##### File Summary <a id="ia-file"></a>
`File Summary` provides a glimpse into the submitted data. This subtab includes
three main components.

1. The `count data` section will detail the first and last five rows of the
   count data and also include the total number of IDs and samples:

   <img src="../vignettes/img/subqc-10.png"/>

2. The second portion includes an overview of the condition data this is
   simply an `R` console-based output of this submitted CSV file:

   <img src="../vignettes/img/subqc-11.png"/>

3. Finally, pre- and post-filtered gene ID counts gives a "before-and-after"
   count of the number of IDs that were filtered using the filter cutoff
   parameter under `2. Data Processing`:

   <img src="../vignettes/img/subqc-12.png"/>

<br/>

##### Count Summary <a id="ia-count"></a>
`Count summary` provides three interactive, downloadable visualizations based
on the file for each sample.

1. Box-and-Whisker Plot for transformed read counts by sample with
   interactivity showing the quartiles, maximum, and minimum count. *With the
   example data, it appears that the box-and-whisker plot for each sample is
   similar to the other samples. If one sample had a plot varying greatly
   from the others, it would indicate some required investigation into that
   specific sample in terms of the number of raw reads provided and proportion
   of reads aligned.*

   <img src="../vignettes/img/subqc-13.png"/>

2. Count Data Distribution plot showing the frequency of transformed count
   data by sample with interactivity displaying the value and frequency.
   Additionally, double-clicking a sample ID in the legend isolates just that
   sample's histogram. Additional sample IDs can be select for more specific
   comparisons also. *With the example data, the histograms appear similar for
   each sample, indicating no significant derivation. Similar to the
   box-and-whisker plots, a single sample varying greatly from other samples
   may indicate a required investigation into the raw read counts or read
   alignment statistics.*

   <img src="../vignettes/img/subqc-14.png"/>
   <img src="../vignettes/img/subqc-15.png"/>


3. Total Reads displays a histogram of total reads by sample with
   interactivity for displaying actual total read counts for each sample.
   Double-clicking on a sample ID in the legend isolates that sample's read
   count histogram and allows for selecting of specific adjacent sample IDs
   for comparison. Total reads counts for individual samples that vary
   greatly from the other total read counts may indicate some issues with data
   preparation (sequencing) or read alignment.  *Here, sample H26313 has a much
   lower total reads count than the other samples. This may be reflected in
   further comparative analyses.*

   <img src="../vignettes/img/subqc-16.png"/>

<br/>

#### Discovery-Driven Analyses <a id="ia-prelim-analysis"></a>
1.  After examining your results on the `Submit and QC` tab, you may proceed to
    the `Discovery-Driven Analyses` tab:

   <img src="../vignettes/img/prelim-01.png"/>

2.  The `Discovery-Driven Analyses` tab will be populated with several subtabs,
    similar to the `Submit and QC` tab. The first subtab you will see is the
    `Correlation` tab. This tab provides correlation analysis of the expression
    matrix by sample:

   <img src="../vignettes/img/prelim-02.png"/>

3.  Under the `Correlation` subtab you will see several visualizations.
    `Interactive Correlation Analysis` displays a heatmap of the correlation
    between samples with interactivity showing the actual correlation between
    the two intersecting samples. *The example data shows most sample-sample
    correlations of 0.95 or larger, indicating relatively high correlation. The
    darker cells here signify less similar samples, which may yield more
    interesting differential expression results. This graph can indicate
    comparisons of interest in future analyses. In most cases, the high number
    of gene IDs with similar or identical expression estimates will cause
    correlations to be large, even for dissimilar genetic expression profiles:*

   <img src="../vignettes/img/prelim-03.png"/>

4.  Clicking on a cell will provided a scatterplot of the two intersecting
    samples with each data point representing one gene ID. A scatterplot with
    points falling more closely along the diagonal indicates samples with more
    similar genetic expression profiles. *This scatterplot shows a clear trend
    of data points (gene IDs) falling along or close to the diagonal. That
    means these two samples have very similar genetic expression profiles. Data
    points that fall far from the diagonal line represent genes with dissimilar
    expression levels between the select samples:*

   <img src="../vignettes/img/prelim-04.png"/>

5.  The Sample Distance Matrix provides a heatmap of the Euclidean distance
    between the gene expression vectors for each sample pair. The larger
    distance (darker red color) indicates samples with the most dissimilar
    genetic expression profiles. This matrix also includes a clustering of the
    samples based on the vectorized expression profiles. *With the example data,
    two distinct clusters can be observed through the first branching of the
    dendrogram. Additionally, as with the correlation heatmap, specific cells
    with a darker color indicate a more dissimilar pair of samples based on
    genetic expression:*

   <img src="../vignettes/img/prelim-05.png"/>

6.  The next visualization will be under the `PCA` subtab. This subtab
    provides Principal Component Analysis (PCA) for the expression estimation
    matrix:

   <img src="../vignettes/img/prelim-06.png"/>

7.  This analysis has the option of selecting a factor of interest. *With the
    example data, selecting “Rootstock” as the factor of interest provides a
    visualization of the first two components for each sample. In this
    application, PCA is a linear transformation of the gene expression levels,
    with the first component representing the transformed dimension with the
    most variability, and each subsequent component decreasing in variability.
    This analysis has the potential to isolate samples based on expression
    levels. Here, there does not appear to be any specific rootstock that
    separates from the others. If there were, it could help develop directions
    for further analysis. The axis labels indicate the first principal
    component accounts for 37% of the variance between samples, whereas the
    second principal component accounts for 7%:*

   <img src="../vignettes/img/prelim-07.png"/>

8.  The next visualization will be under the `MDS` subtab. This subtab
    provides multi dimensional scaling for the expression estimation matrix.
    This is similar to PCA except is develops components through non-linear
    transformations of the gene expressions:

   <img src="../vignettes/img/prelim-08.png"/>

9.  *Looking back at our sample data, we observe similar results to that of the
    PCA results, with similar potential interpretations if any sample or groups
    of samples were to differentiate from the others:*

   <img src="../vignettes/img/prelim-09.png"/>

10. The next visualization will be under the `Heatmap` subtab. This subtab
    provides an interactive heatmap with rows representing the gene IDs and
    columns representing sample IDs:

   <img src="../vignettes/img/prelim-10.png"/>

11. This heatmap requires an ID cutoff to select the indicated number of gene
    IDs with the highest mean expression values for display. Selecting a cell
    displays a plot showing the total read counts for that specific gene ID by
    the selected factor. *With the example data, the 20 most variable gene IDs
    are displayed. The yellow color indicates gene IDs with a higher
    expression level for that sample, and the darker blue color represents a
    low expression level for that sample. Selecting ID: rna25007 shows the
    read counts for that ID by rootstock factor. This shows the “OWN”
    rootstock seems to have a higher expression level for that ID, with the
    exception of one sample:*

   <img src="../vignettes/img/prelim-11.png"/>
   <img src="../vignettes/img/prelim-12.png"/>

12. On the next subtab (1), `Biclustering` performs a biclustering analysis
    using one of the selected biclustering tools (3) with a maximum bicluster size
    of the indicated cutoff value (2):

   <img src="../vignettes/img/prelim-13.png"/>

13. Launching the analysis results in display of the first bicluster.
    Alternative clusters can be selected from the dropdown menu and the IDs
    and plot for each bicluster can be downloaded using the button below the
    visualization. *With the sample data, the biclusters can help select the
    samples under which certain gene IDs are similarly expressed. Since gene
    expression levels can vary greatly over all samples and conditions, a
    biclustering approach can isolate similar expression patterns on a level
    where traditional clustering may miss. The first cluster for the example
    data shows that samples B20715, D21515, and H12915 are expressed similarly
    under the isolated gene IDs. Interpretations can be made similarly for
    each subsequent bicluster:*

   <img src="../vignettes/img/prelim-14.png"/>

<br/>

#### DGE Analysis <a id="ia-dge-analysis"></a>
1.  After uploading the user data or selecting the example data, users can go
    directly to the `DGE Analysis` tab, preceding the `Discovery-Driven Analyses`
    functions:

   <img src="../vignettes/img/dge-01.png"/>

2.  Once you are on this tab, you will be greeted with several options.
    `1. Experimental Setup` will allow users to select an experimental design
    for their DGE analysis. Options are Two-group comparison, Multiple factor
    comparisons, Classical interaction design, Additive models (paired or
    blocking designs), Main effects, Main effects with grouping factors, or
    Custom design:

   <img src="../vignettes/img/dge-02.png"/>

<br/>

##### An overview of experimental designs <a id="ia-over-exp"></a>
1.  The `Two-group comparisons` options is the traditional approach for DGE and
    compares two factors levels for the selected factor type. *With the example
    data, selecting “Two group comparisons” for the experimental design and
    “Rootstock” for the factor allows for specific pairwise comparisons of
    Rootstock factor levels. Here, we can select specific comparisons of
    interest from the permutations of all pairwise comparisons. Selecting all
    comparison options will provide inverse duplications, so specific
    selections may be needed. Below, all unique pairwise combinations are
    selected. The linear model is also displayed for users interested in the
    model used for comparisons:*

   <img src="../vignettes/img/dge-03.png"/>

2.  The `Multiple factor comparisons` design has users select two factor levels
    and performs all crosswise comparisons for the two chosen factor levels.
    *With the example data, the Multiple factor comparisons design with
    Rootstock and Block selected as the two factors provides optional
    comparisons for each rootstock separated by block. In this situation, as
    with the other designs, the user selects which comparisons are of interest.
    Selecting C3309_B_VS_C3309_E allows for a comparisons of gene expression
    levels for the same rootstock in two different blocks. This provides
    insight into the locations and possibly time (due to time requirements for
    sampling) for this specific rootstock:*

   <img src="../vignettes/img/dge-04.png"/>

3.  The `Classical interaction design` allows the user to select two factor
    levels and one reference level for each of the selected factors. *With the
    experimental data, the Classic interaction design allows for selection of
    two factors and a reference for each. In this case, the Rootstock and Row
    are the chosen factors and OWN and 15 are the selected reference levels for
    comparison. The contrast levels selected upon submission will provide DEGs
    with respect to these two levels:*

   <img src="../vignettes/img/dge-05.png"/>

4.  The `Additive models` design is useful when samples are paired or blocked
    into distinct groupings. This format requires selected factors, one for
    the pairing or blocking factor and the other for the treatment factor.
    Additionally, a reference level for each selected factor is required. *In
    the example data, the Rootstock factor can be considered for grouping and
    the treatment factor is Block. OWN is selected as the reference level for
    the blocking factor and the A15 Block is selected for the treatment
    reference level:*

   <img src="../vignettes/img/dge-06.png"/>

5.  `Main Effects` experimental design allows for testing the significance of a
    factor across multiple factor levels. In this situation, any significant
    deviation from an intercept would result in a significantly differential
    expressed gene. Using this approach is most useful when users want to test
    the significance of a factor that has more than two levels. This design
    requires indication of which factor to test as the main effect.
    Additionally, users must specify a factor reference level for the lfc and
    for corresponding visualizations. The $p$-values, adjusted $p$-values, and
    related differentially expressed genes are not affected by the selected
    reference level. *In the example data, the Rootstock factor is selected,
    which will test whether Rootstock has a significant impact on the genetic
    expression across all Rootstock factor levels. Selection of the C3309
    Rootstock factor level as a reference will provide lfc values and
    visualizations relative to this level, while $p$-values and adjusted
    $p$-values will be calculated with respect to the main effect of the
    chosen factor:*

   <img src="../vignettes/img/dge-07.png"/>

6.  `Main effects with grouping factors` design allows for a more detailed main
    effects testing to be performed. This design tests the indicated main
    effect when the data is subset by a user-indicated factor. Requirements for
    this design are indication of main effect to test, factor for grouping, and
    grouping factor level for which to provide results. As with the main
    effects design, the user also specifies a main effect reference level for
    lfc and visualizations. *In the example data, Rootstock is selected as the
    grouping factor, with C3309 selected as the grouping factor level for which
    the results table will be based on. The Block factor is selected to test
    the main effect, with B13 selected as the main effect reference level. This
    analysis will indicate which genes are differentially expressed when
    account for Block only for the C3309 samples:*

   <img src="../vignettes/img/dge-08.png"/>

7.  Finally, `Custom design` provides advanced users with a method for
    indicating their own design matrix. This method is provided so that users
    can perform more intricate DGE analyses beyond what IRIS-EDA already
    provides as experimental design options. It is only recommended that users
    with advanced knowledge of the modeling process to use this option, as the
    required design matrix may not be intuitive to use appropriately:

   <img src="../vignettes/img/dge-09.png"/>

<br/>

##### Setting up an experiment <a id="ia-set-exp"></a>
1.  Following experimental design selection, the user can then select the
    desired tool for performing differential gene expression analysis. This
    will be found under `2. DGE Parameters`. Options for this purpose are
    `DESeq2`, `edgeR`, and `limma-voom`. If the user selects one of the main
    effects design options, the only options for DGE tools are `DESeq2` and
    `edgeR`:

   <img src="../vignettes/img/dge-10.png"/>

2.  Also under `2. DGE Parameters`, adjusted $p$-value and log-fold change
    cutoffs can be specified to filter the results provided. The default values
    are 0.05 for adjusted $p$-value and 1 for the minimum log-fold change. To
    provide all gene comparisons, select 1 for adjusted $p$-value cutoff and 0
    for the minimum log-fold change:

   <img src="../vignettes/img/dge-11.png"/>

3.  Upon indication of the experimental design, DGE tool and cutoffs, the
    parameters and comparisons for the specified analysis must be provided.
    Such parameters can be found under `3. Experimental Parameters`. Depending
    on the selected experimental design, the required parameters and
    comparisons will vary:

   <img src="../vignettes/img/dge-12.png"/>

4.  Launching the analysis will perform the DGE analysis and provide results.
    This is found under `4. Launch Analysis`. **Note:** Be mindful that this
    process may take more time than previous steps. Adjusting the cutoff
    values after submission will automatically update the generated
    visualizations and results table:

   <img src="../vignettes/img/dge-13.png"/>

5.  If you haven't crashed your computer or our server, the main page will now
    be populated with data in the two subtabs: `Overview` and `Plots`.

6.  The `Overview` subtab (1) will display information based on the number of
    differentially expressed genes (DEGs), including a table of significantly
    up- and down-regulated gene ID counts by selected comparisons (2) and an
    interactive barplot representing this data (3). *With the example data, the
    results show only one gene is differentially expressed in either direction
    from the two-group comparison based on the Rootstock factor. This
    information seems to follow with the previously investigated figures
    showing limited clustering of samples from the PCA and MDS and high
    correlation values between samples:*

   <img src="../vignettes/img/dge-14.png"/>

7.  The `Plots` subtab (1) will provides interactive visualizations based on the
    differential expression analysis as well as the actual results file from
    the previously-selected tool based on the selected comparison. The MA plot (2)
    shows the transformed log-fold change compared to the transformed base
    mean for each gene ID. Specific points can be selected, highlighting both
    the point and the corresponding row in the results table found below the
    figure. Additionally, the gene ID can be selected from the results table,
    which will highlight the location of that gene ID on the figure:

  <img src="../vignettes/img/dge-15.png"/>

8.  The Volcano plot shows a comparison of the transformed $p$-value compared
    with the transformed log-fold change. As with the MA plot, the Volcano plot
    is interactively connected with the results table found below the figure:

  <img src="../vignettes/img/dge-16.png"/>

9. Finally, the results table is generated in coordination with the above plots
   and displays the output file of the selected tool’s differential gene
   expression analysis. This table can be sorted increasingly or decreasingly
   by any of column. The integrated search feature also allows for specific
   gene IDs to be found in the table. This results file can be exported in a
   filtered or unfiltered format. *With the example data, sorting by adjusted
   $p$-value (padj) shows that there are no differentially expressed genes for
   this specific comparison (C3309 vs. OWN):*

  <img src="../vignettes/img/dge-17.png"/>

<br/>

<br/>

### GEO Usage <a id="geo-main"></a>

#### Overview <a id="geo-over"></a>
GEO (**G**ene **E**xpression **O**mnibus) is a public data repository that
accepts array- and high throughput sequence-based data. To submit data to GEO,
you will need three components:

  * Metadata spreadsheet
  * Processed data files
  * Raw data files

The `GEO` section of this web server will aid in the production of two
of the prior components: **metadata and processed data generation.**
This works by having the user fill out a dynamic questionnaire
which follows the template foundation located on NCBI's servers shown
[here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html). While
the user fills out this form, IRIS-EDA will automatically parse your
raw count data into seperate text files based on the number of samples
found in your raw count matrix. This can be downloaded via ZIP file.
Additionally, the information provided by the column data matrix will
be used to populate several components of the metadata, including
sample information, characteristics, and
[MD5 checksum generation](https://en.wikipedia.org/wiki/MD5).

To get started, you *must* load both the raw count data and sample
information on the `Submit and QC` tab. If you do not and head
directly to the `GEO` tab, you will be greeted with a page that looks
like this:

<br/>

<img src="../vignettes/img/geo-series-01.png"/>

<br/>

#### Series <a id="geo-series"></a>
After you have successfully submitted your sample data and sample
information, the first section that you will need to fill out is the
**Series** section. This section describes the overall experiment.
Specifically, there are six major components:

  1.  **Title**: a unique title (*less than 255 characters*) that
      describes the overall study.

  2.  **Summary**: a *thorough* description of the goals and objectives
      of this study. The abstract from the associated publication may
      be suitabe. You may provide as much text here as necessary.

  3.  **Overall design**: Indicate how many samples are analyzed, if
      replicates are included, are there control and/or reference
      Samples, etc.

  4.  **Contributor(s)**: Who contributed to this study? You must
      provide the information in the following format
      (`Firstname,Initial,Lastname`). Examples of this format:

        * "Brandon,T,Monier"
        * "Jane,Doe"

      This component is dynamic in this framework. If there was more
      than one contributor, click the `Add Contributor` button to add
      as many individuals who aided in this endeavor.
  5.  **Supplementary file**: If you want to submit your raw count
      matrix provided earlier, include the file name here (**optional**).

  6.  **SRA center name code**: Only enter a value if your institute
      already has a `Center_Name` code registered with SRA. Otherwise,
      leave empty (**optional**).

<br/>

#### Samples <a id="geo-samples"></a>
This section lists and describes each of the biological Samples
under investgation, as well as any protocols that are specific
to individual Samples. Additional "processed data files" or "raw files"
may be included.

Additionally, this section is dynamically generated. The tabs that will
populate this section will be based on the number of samples you have
submitted for analysis. For example, the "small example data set"
has seven samples, therefore the number of tabs that would have to
be filled out would be seven:

<br/>

<img src="../vignettes/img/geo-samples-01.png"/>

<br/>

**It is mandatory that you fill out all samples.** In order for the
spreadsheet to properly be formatted, all of these sample tabs need to
be filled out. Each sample has the following criteria:

  1.  **Sample name**: This will automatically be filled out for you.
      This provides an unique identifier for each sample that will
      not affect any downstream naming schemes.

  2.  **Title**: A unique title that describes the sample.

  3.  **Source name**: Briefly identfiy the biological material.
      (e.g. vastus lateralis muscle)

  4.  **Organism**: Identify the organism(s) from which the
      sequences were derived.

  5.  **Molecule**: The type of molecule that was extracted from the
      biological material. Choose from one of the preselected items
      from the dropdown input.

  6.  **Description**: Additional information not provided in the
      other fields, or paste in broad descriptions that cannot be
      easily dissected into the other fields.

  7.  **Processed data file**: The name of the file containing the
      processed data. **Note (1)**: This step has partially been
      completed for you. This will be the individual sample text
      file of raw counts. If you have more than one processed data
      file that you would like to contribute to GEO, click the
      `Add Proc. Data File` button to add as many entries you need.

  8.  **Raw file**: The name of the file containing the raw data.
      Similar to the processed data, additional raw data can be
      provided by clicking the `Add Raw Data File` button.

**Note (2)**: If you have taken a look at the metadata template before,
you will also notice that this section contains a characteristic column.
This column is automatically generated for you based on your sample
information provided in conjuntion with your raw count matrix. For example,
the "small example data set" contains the following columns of conditions:

  * condition
  * type

Therefore, **two** characteristics columns will be generated in the
spreadsheet. Each of the rows in these columns will based on the factor
levels for each condition.

<br/>

#### Protocols <a id="geo-protocols"></a>
This section provides the public with information about what protocols
were used to conduct this experiment. This section has the following criteria:

  1.  **Growth protocol**: Describe the conditions that were used to
      grow or maintain organisms or cells prior to extract
      preparation (**optional**).

  2.  **Treatment protocol**: Describe the treatments applied to the
      biological material prior to extract preparation.

  3.  **Extract prototocol**: Describe the protocols used to extract and
      prepare the material to be sequenced.

  4.  **Library construction protocol**: Describe the library construction
      protocol.

  5.  **Library strategy**: A Short Read Archive-specific field that
      describes the sequencing technique for this library. Choose from one
      of the preselected items from the dropdown input.

<br/>

#### Data processing pipeline <a id="geo-pipe"></a>
This section includes steps for various data processing techniques.
Base-calling, alignment, filtering, peak-calling, generation of normalized,
abundance measurements, etc.

  1.  **Data processing step**: Provide details of how processed data were
      generated. This can include procedures described in the paragraph
      above. Similar to prior steps, this section is dynamic. To add more
      steps, click the `Add Data Proc. Step` butto.

  2.  **Genome build**: UCSC or NCBI genome build number (e.g. hg18,
      mm9, etc.) or reference sequence used for read alignment.

  3.  **Processed data files format and content**: For each processed
      data file type, provide a description of the format and content. 

<br/>

#### Processed data files <a id="geo-proc-data"></a>
For each file that you provided additional processed data in the 
`Samples` section, you will need to provide additional information about 
these files below:

  1.  **File type**: The type of processed file. Examples include, peak,
      wig, bed, gff, bigWig, etc.

  2.  **MD5 file checksum**: MD5 checksum of the file. This helps GEO
      verify that the file transfer was complete and didn't corrupt your 
      file.

<br/>

#### Raw files <a id="geo-raw-data"></a>
Similar to the `Processed data files` section above, for each file listed
in the "raw files" of the `Samples` section, provide additional information
about these files below:

  1.  **File name**: The name of the file provided in the `Samples` section.

  2.  **File type**: The type of processed file. Examples include, peak,
      wig, bed, gff, bigWig, etc.

  3.  **MD5 file checksum**: MD5 checksum of the file. This helps GEO
      verify that the file transfer was complete and didn't corrupt your 
      file.

  4.  **Instrument model**: Include the instrument make and model used to
      sequence the samples. Choose from one of the preselected items from 
      the dropdown input.

  5.  **Read length**: The number of bases expected in each raw sequence. 
      If you are using variable read length technology (e.g. 454, PacBio,
      Ion Torrent), put `0` for variable read length.

  6.  **Single or paired-end**: Choose either if your samples were sequenced
      with single reads or paired-end reads. If you used SOLiD technology,
      choose `SOLiD` in the options. **This is critical for the final 
      section**.

<br/>

#### Paired-end experiments / SOLiD data <a id="geo-pair"></a>
This section will only need to be filled out if your raw files used
paired-end or SOLiD-based sequencing technologies. For paired-end experiments,
list the 2 associated raw files, along with some additional information:

  1.  **File name 1**: Paired-end read 1.

  2.  **File name 2**: Paired-end read 2.

  3.  **Average instert size**: Average size of the insert for paired-end
      reads.

  4.  **Standard deviation**: Standard deviation of insert size. This is
      typically around 10% of the insert size stated.

**Note (1)**: The previously mentioned components are for standard 
paired-end read sets only. If your raw reads used SOLiD technology, this
section will may appear void of any information (see example):

<br/>

<img src="../vignettes/img/geo-pair-01.png"/>

<br/>

Instead of 2 reads, you will need to provide 4 (due to the basis of this 
technology) In the final section before the downloads portion. Therfore, 
instead of 2 "File names", You will be prompted to fill out 4; each 
corresponding to their respective read.

<br/>

#### Downloads <a id="geo-downloads"></a>
After you have submitted all of your information, you can download two items:
  
  * An Excel spreadsheet of the metadata;
  * A ZIP file of individual sample reads as text files.

<br/>

<br/>

### Sing-cell RNA-seq <a id="scrna-main"></a>

#### Overview <a id="scrna-over"></a>
Upon sample submission to IRIS-EDA, you have two options: 
(1) Perform analyses on a normal RNA-seq experiment or 
(2) run single-cell RNA-seq analysis. 

If you choose the latter option, you will need to provide a file 
containing ID lengths to help reduce variability. This is crucial since
filtration of low TPM (transcripts per million) reads can lead to better 
results. In order to calculate TPM, ID length approximations are needed
in conjunction with raw count data. To submit this information to IRIS-EDA, 
**you will need to provide ID data that follows the proceeding criteria**:

```
id        id_len_kbp
gene001          1.4
gene002          5.6
gene003          0.8
...              ...
```

**Note**:
This file must be in a comma seperated value (CSV) format, where column 1
is the ID (e.g. gene, transcript, etc.) and column 2 is the approximation
of length in **kilobase pairs**.

<br/>

#### A method to obtain this data <a id="scrna-methods-01"></a>
There are many ways to obtain these values. A common procedure would be to
parse the respective general fearture format version 3 (GFF3) file and 
determine the length through calculating the difference between the 
start and end locations. 

If you need help parsing this information, a primitive `R` function has been
made, which you can find 
[here](https://gist.github.com/btmonier/409856fe22280603ca19345fdd4e5e22).
To run this function, you will need to install and load 
3 packages into your R library:

  * `ape`
  * `stringr`
  * `dplyr`

All of these packages can be found on the 
[CRAN repository](https://cran.r-project.org/), or by using the
`install.packages()` function provided in the base package of `R`. 

This function has four parameters:

  * `gff`: the location to the GFF3 file on your computer

  * `cts`: your raw count matrix that you will provide to IRIS-EDA. **Note**:
    this object must be a `matrix` type and have the same format that is
    found in the first section of this walkthrough (*A note about input data*).

  * `type`: the GFF3 type column identifier (e.g. gene, transcript, exon,
    CDS, etc.)

  * `attribute`: the GFF3 attribute (e.g. `ID`, `gene_id`, 
    `transcript_id`, `gene_name`, or `transcript_name`)

Depending on the size of this file, this may take some time. This function
will return a `tibble` data frame, so make sure you assign it to an object
before you use `write.csv()`. After you perform this task, the file can
successfully be submitted to IRIS-EDA. 

<br/>

#### TPM filtration <a id="scrna-tpm-filter"></a>
Once you have submitted the data, you will notice that the `Filter cutoff`
changes from `count data row sums` to `TPM`:

<br/>

<img src="../vignettes/img/geo-tpm-01.png"/>

<br/> 

The default is set to a value of `1`, however, this can be changed at the
user's discretion. **Note (1)**: this value corresponds to which rowsums
(i.e. ID sums) will be filtered out if they have a value that is less than
the user parameter.

<br/>

<br/>
