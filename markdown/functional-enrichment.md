<br/>
Since there are many robust and capable servers for functional enrichment
analysis, we feel it is better to use a previously tool/web server than to 
try to reproduce what is already available.

Examples of these services include the following:
  * [DAVID](https://david.ncifcrf.gov/)
  * [UCSC Genome Browser](https://genome.ucsc.edu/)
  * [UniProt](https://www.uniprot.org/)

For this example, we will look at how you can extract differentially 
expressed gene data for IRIS-EDA and how you can use DAVID for functional
enrichment.

1. Once you have run **DGE Analysis** in the `DGE Analysis` tab, simply
   download the differentially expressed genes using the 
   `Download Filtered Data` button at the bottom of the page on the 
   prior tab (`Plots` tab).
   
2. This will download a set of *significantly* differentially expressed genes 
   or IDs that can be opened up in R, Excel, LibreOffice Calc, or any other
   spreadsheet viewing software.
   
3. Once the spreadsheet is open, simply copy the IDs from the first column and
   open up DAVID Bioinformatics Resources (https://david.ncifcrf.gov/).
   
4. Click on the `Gene Functional Classification` button on the **left-hand** 
   side of the page. *This will be under the Shortcut to DAVID Tools*.
   
5. On this page, click on the `Upload` button on the **left-hand** side of the
   page. In *Step 1: Enter Gene List*, paste your gene IDs into this text
   box. *If you have pasted this list into a text file, you can upload 
   that instead.*
   
6. For *Step 2: Select Identifier*, choose the ID source for streamlined 
   analysis. For example, the small example data set found in the tutorial is
   derived from the `FLYBASE_GENE_ID` database, therefore, I would select
   that from the list. Your data, of course, may differ. If you don't know
   the origin of the gene IDs,there is a `Not Sure` option at the bottom of 
   the list.
   
7. For *Step 3: List Type*, choose `Gene List` and finally, submit the list
   using the `Submit List` button in *Step 4: Submit List*. 
   
8. Once submitted, you can view which genes are enriched, clusters, etc.
   For more information about downstream analyses, please click 
   [here](https://david.ncifcrf.gov/helps/functional_classification.html).
