#---------------------------------------------------------------------
# Title:         IRIS - Server Script
# Author:        Brandon Monier
# Created:       2018-01-26 at 11:32:02
# Last Modified: 2018-10-29 at 11:01:10
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# CONTENTS
#---------------------------------------------------------------------
#
#   IRIS server function organization
#
#       01 DATA...... data loading
#       02 QC........ quality control for preliminary analysis
#       03 DGE....... differential gene expression analysis
#       04 HEAT...... heatmap plots generation
#       05 BIC....... biclustering algorithms and visualization
#       06 COR....... correlation analysis and visualization
#       07 SESSION... generate session info
#       08 GEO....... gene expression omnibus metadata generation
#       09 TABS...... tab links using `updateNavbarPage`
#
#---------------------------------------------------------------------
#
#   IRIS application layout
#
#       IRIS
#       |-  Welcome
#       |-  Submit and QC
#           |-  File Summary
#           |-  Count Summary
#       |-  Discovery-Driven Analysis
#           |-  Correlation
#           |-  PCA
#           |-  MDS
#           |-  tSNE
#           |-  Heatmap
#           |-  Biclustering
#           |-  Clustering
#       |-  DGE Analysis
#           |-  Overview
#           |-  Plots
#       |-  GEO
#           |- Dynamic sample generation (see "GEO" layout)
#       |-  More
#           |-  Tutorial
#           |-  FAQ
#           |-  About Us
#           |-  Session Info
#
#---------------------------------------------------------------------
#
#   IRIS gene expression omnibus (GEO) layout
#
#       IRIS
#       |-  GEO
#           |-  Series
#               |-  Title
#               |-  Summary
#               |-  Overall design
#               |-  Contributor (DYNAMIC)
#               |-  Supplementary file
#               |-  SRA center name code
#           |-  Samples
#               |-  Sample name
#               |-  Title
#               |-  Source name
#               |-  Organism
#               |-  Chacteristics (BASED ON METADATA)
#               |-  Molecule
#               |-  Description
#               |-  Processed data file (DYNAMIC)
#               |-  Raw file (DYNAMIC)
#           |-  Protocols
#               |-  Growth protocol
#               |-  Treatment protocol
#               |-  Extract protocol
#               |-  Library construction protocol
#               |-  Library strategy
#           |-  Data processing pipeline
#               |-  Data processing step (DYNAMIC)
#               |-  Genome build
#               |-  Processed data files format and content
#           |-  Processed data files
#               |-  File name
#               |-  File type
#               |-  File checksum
#           |-  Raw files
#               |-  File name
#               |-  File type
#               |-  File checksum
#               |-  Instrument model
#               |-  Read length
#               |-  Single or paired-end
#           |-  Paired-end experiments (IF RAW FILE -> PAIRED END)
#               |-  File name 1
#               |-  File name 2
#               |-  Average insert size
#               |-  Standard deviation
#           |-  SOLiD experiments (IF RAW FILE -> SOLiD)
#               |-  File name 1
#               |-  File name 2
#               |-  File name 3
#               |-  File name 4
#               |-  Average insert size
#               |-  Standard deviation
#           |-  Download options
#               |-  Metadata (.xlsx)
#               |-  Processed data (.zip)
#               |-  Both meta- and processed data (.zip)
#
#---------------------------------------------------------------------


# Change file upload size to 30 MB and sanitize errors
options(
    shiny.maxRequestSize = 30 * 1024^2,
    shiny.sanitize.errors = TRUE
)



# Server function
irisServer <- function(input, output, session) {

    ###################################################################
    ###################################################################
    ### SECTION 01 - DATA LOAD (DATA)
    ###################################################################
    ###################################################################

    ## DATA - Data load option - count data
    output$file1 <- renderUI({
        if (input$examplechoice == "no" | input$examplechoice == "scrna" | input$examplechoice == "scrna10x") {
            fileInput(
                inputId = "file1",
                label = "Submit count data (CSV)",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            )
        } else {
            return()
        }
    })

    ## DATA - Data load option - metadata
    output$file2 <- renderUI({
        if (input$examplechoice == "no" | input$examplechoice == "scrna" | input$examplechoice == "scrna10x") {
            fileInput(
                inputId = "file2",
                label = "Submit metadata (CSV)",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            )
        } else {
            return()
        }
    })

    ## DATA - Data load option - gene length
    output$file3 <- renderUI({
        if (input$examplechoice == "scrna") {
            fileInput(
                inputId = "file3",
                label = "Submit ID lengths (CSV)",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            )
        } else {
            return()
        }
    })

    ## DATA - Filter choice
    output$filter_choice <- renderUI({
        if (input$examplechoice == "yes3" | input$examplechoice == "scrna") {
            textInput(
                inputId = "prefilt",
                label = withMathJax(
                    "Filter cutoff (TPM \\( < n\\))"
                ),
                value = as.numeric(1)
            )
        } else if (input$examplechoice == "scrna10x") {
            textInput(
                inputId = "prefilt",
                label = withMathJax(
                    "Filter cutoff (TPM \\( < n\\))"
                ),
                value = as.numeric(0)
            )
        } else {
            textInput(
                inputId = "prefilt",
                label = withMathJax(
                    "Filter cutoff (count data row sums \\( < n\\))"
                ),
                value = as.numeric(10)
            )
        }
    })





    ###################################################################
    ###################################################################
    ### SECTION 02 - QUALITY CONTROL (QC)
    ###################################################################
    ###################################################################

    ## QC - reactive - load and add data to DESeqDataSet class
    ddsout <- eventReactive(input$goqc, {
        if (input$examplechoice == "yes1") { # Small example
            cts <- as.matrix(
                read.csv(
                    "./data/count-data-small.csv",
                    header=TRUE,
                    row.names = 1
                )
            )
            coldata <- read.csv(
                "./data/col-data-small.csv",
                header = TRUE,
                row.names = 1
            )
            cols <- colnames(coldata)
            coldata[cols] <- lapply(coldata[cols], factor)
            cts <- cts[, rownames(coldata)]
            n <- input$prefilt
            n <- as.numeric(n)
            cts.filt <- cts[rowSums(cts) > n, ]
            dds <- DESeqDataSetFromMatrix(
                countData = cts.filt,
                colData = coldata,
                design = ~ 1
            )
        } else if (input$examplechoice == "yes2") { # Big example
            cts <- as.matrix(
                read.csv(
                    "./data/count-data-big.csv",
                    header = TRUE,
                    row.names = 1
                )
            )
            coldata <- read.csv(
                "./data/col-data-big.csv",
                header = TRUE,
                row.names = 1
            )
            cols <- colnames(coldata)
            coldata[cols] <- lapply(coldata[cols], factor)
            cts <- cts[, rownames(coldata)]
            n <- input$prefilt
            n <- as.numeric(n)
            cts.filt <- cts[rowSums(cts) > n, ]
            dds <- DESeqDataSetFromMatrix(
                countData = cts.filt,
                colData = coldata,
                design = ~ 1
            )
        } else if (input$examplechoice == "yes3") { # scRNA example
            cts <- read.csv(
                "./data/count-data-scrna.csv",
                header = TRUE,
                row.names = 1
            )
            coldata <- read.csv(
                "./data/col-data-scrna.csv",
                header = TRUE,
                row.names = 1
            )
            gen.len <- read.csv(
                "./data/gen-len-scrna.csv",
                header = TRUE,
                row.names = 1
            )

            cols <- colnames(coldata)
            coldata[cols] <- lapply(coldata[cols], factor)

            cts.nam <- row.names(cts)
            gen.len <- gen.len[order(
                match(row.names(gen.len), row.names(cts))), , drop = FALSE]

            # Normalize for gene length (TPM)
            x <- cts / as.matrix(gen.len)

            # Normalize for sequencing depth (TPM)
            tpm.mat <- t( t(x) * 1e6 / colSums(x) )

            # Filter raw counts based TPM values
            n <- input$prefilt
            n <- as.numeric(n)
            tpm.mat <- cts[rowSums(cts) > n, ]
            cts.filt <- cts[rownames(tpm.mat), ]

            # Assign to DESeqDataSet object
            dds <- DESeqDataSetFromMatrix(
                countData = cts.filt,
                colData = coldata,
                design = ~ 1
            )

        } else if (input$examplechoice == "no" | input$examplechoice == "scrna10x") { # user input (normal)
            cts <- input$file1
            coldata <- input$file2
            cts <- read.csv(
                cts$datapath,
                header = TRUE,
                row.names = NULL
            )
            #set the first column name to 'gene_id' if missing
            if(colnames(cts)[1] == ""){
              colnames(cts)[1] == "gene_id"
            }
            # if there are duplicated row names, only keep the first
            if (length(which(duplicated.default(cts[,1]))) > 0) {
              cts <- cts[-which(duplicated.default(cts[,1]) == T),]
            }
            rownames(cts) <- cts[,1]
            cts <- as.matrix(cts[,-1])
            coldata <- read.csv(
                coldata$datapath,
                header = TRUE,
                row.names = 1
            )
            cols <- colnames(coldata)
            coldata[cols] <- lapply(coldata[cols], factor)
            cts <- cts[, rownames(coldata)]
            n <- input$prefilt
            n <- as.numeric(n)
            cts.filt <- cts[rowSums(cts) > n, ]
            dds <- DESeqDataSetFromMatrix(
                countData = cts.filt,
                colData = coldata,
                design = ~ 1
            )
        } else if (input$examplechoice == "scrna") { # user input (scRNA)
            cts <- input$file1
            coldata <- input$file2
            gen_len <- input$file3
            cts <- read.csv(
                cts$datapath,
                header = TRUE,
                row.names = NULL
            )
            #set the first column name to 'gene_id' if missing
            if(colnames(cts)[1] == ""){
              colnames(cts)[1] == "gene_id"
            }
            # if there are duplicated row names, only keep the first
            if (length(which(duplicated.default(cts[,1]))) > 0) {
              cts <- cts[-which(duplicated.default(cts[,1]) == T),]
            }
            rownames(cts) <- cts[,1]
            cts <- as.matrix(cts[,-1])
            coldata <- read.csv(
                coldata$datapath,
                header = TRUE,
                row.names = 1
            )
            gen_len <- read.csv(
                gen_len$datapath,
                header = TRUE,
                row.names = 1
            )

            cts.nam <- row.names(cts)
            gen_len <- gen_len[order(
                match(row.names(gen_len), row.names(cts))), , drop = FALSE]

            # Normalize for gene length (TPM)
            x <- cts / as.matrix(gen_len)

            # Normalize for sequencing depth (TPM)
            tpm.mat <- t( t(x) * 1e6 / colSums(x) )

            # Filter raw counts based TPM values
            n <- input$prefilt
            n <- as.numeric(n)
            tpm.mat <- cts[rowSums(cts) > n, ]
            cts.filt <- cts[rownames(tpm.mat), ]

            # Assign to DESeqDataSet object
            dds <- DESeqDataSetFromMatrix(
                countData = cts.filt,
                colData = coldata,
                design = ~ 1
            )
        }

        # Return objects for downstream analysis
        return(list(dds, coldata, cts.filt, cts))
    })

    ## QC - reactive - data transformation
    ddstran <- eventReactive(input$goqc, {
        dds <- ddsout()[[1]]
        if (input$transform == "log") {
            tmp <- normTransform(dds)
            lab <- "log<sub>2</sub>(counts + 1)"
        } else if (input$transform == "rlog") {
            tmp <- rlog(dds)
            lab <- "rlog(counts)"
        } else if (input$transform == "vst") {
            tmp <- vst(dds)
            lab <- "vst(counts)"
        } else if (input$transform == "raw") {
            tmp <- dds
            lab <- "Raw counts"
        }
        return(list(tmp, lab))
    })

    ## QC - header (2) - file summary (count data)
    output$filesummarycts <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Count data (first 5 rows and last 5 rows)")
        }
    })

    ## QC - header (2) - file summary (col data)
    output$filesummarycoldata <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            h4("Sample metadata")
        }
    })

    ## QC - header (2) - ID counts (pre)
    output$headcountpre <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            h4("Number of IDs (pre-filtration)")
        }
    })

    ## QC - header (2) - ID counts (post)
    output$headcountpost <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            h4("Number of IDs (post-filtration)")
        }
    })

    ## QC - verbatim console - count data head
    output$fileoutputcts <- renderPrint(width = 400, {
        cts <- ddsout()[[4]]
        trubble(cts)
    })

    ## QC - verbatim console - metadata
    output$fileoutputcoldata <- renderPrint({
        coldata <- ddsout()[[2]]
        coldata
    })

    ## QC - verbatim console - gene no. (pre-filter)
    output$fileoutputcountpre <- renderPrint({
        cts.pre <- ddsout()[[4]]
        nrow(cts.pre)
    })

    ## QC - verbatim console - gene no. (pre-filter)
    output$fileoutputcountpost <- renderPrint({
        cts.post <- ddsout()[[3]]
        nrow(cts.post)
    })

    ## QC - header (2) - boxplot
    output$countbox <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Count data distributions - box and whisker")
        }
    })

    ## QC - visualize - boxplot
    output$boxplot <- renderPlotly({
        withProgress(message = "Compiling boxplots...", value = 0, {
            incProgress()
            input$goqc
            tmp <- ddstran()[[1]]
            tmp <- assay(tmp)
            lab <- ddstran()[[2]]
            isolate({
                box <- as.data.frame(tmp)
                box <- tidyr::gather(box)
                plot_ly(
                    box,
                    type = "box",
                    y = ~value,
                    color = ~key
                ) %>%
                layout(
                    margin = list(b = 90),
                    xaxis = list(title = "", tickangle = -45),
                    yaxis = list(title = lab)
                )
            })
        })
    })

    ## QC - DOWNLOAD BUTTON - Boxplot (pdf)
    output$dlqcboxplotpdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcboxplotpdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - DOWNLOAD PLOT - Boxplot (pdf)
    output$dlqcboxplotpdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-boxplot.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6.5, onefile = FALSE) # open the pdf device
            qcBoxPlot(
                tmp = ddstran()[[1]],
                lab = ddstran()[[2]],
                tran = input$transform
            )
            dev.off()
        }
    )

    ## QC - DOWNLOAD BUTTON - Boxplot (png)
    output$dlqcboxplotpng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcboxplotpngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - DOWNLOAD PLOT - Boxplot (png)
    output$dlqcboxplotpngimg <- downloadHandler(
        filename =  function() {
            paste("qc-boxplot.png")
        },
        content = function(file) {
            png(file, width = 900, height = 750)
            qcBoxPlot(
                tmp = ddstran()[[1]],
                lab = ddstran()[[2]],
                tran = input$transform
            )
            dev.off()
        }
    )


    ## QC - header (2) - histogram
    output$counthist <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            h4("Count data distributions - histogram")
        }
    })

    ## QC - visualize - histogram
    output$hist <- renderPlotly({
        withProgress(message = "Compiling histograms...", value = 0, {
            incProgress()
            input$goqc
            tmp <- ddstran()[[1]]
            tmp <- assay(tmp)
            lab <- ddstran()[[2]]
            isolate({
                hist <- as.data.frame(tmp)
                hist <- tidyr::gather(hist)
                plot_ly(
                    data = hist,
                    type = "histogram",
                    x = ~value,
                    color = ~key
                ) %>%
                layout(
                    xaxis = list(title = lab),
                    yaxis = list(title = "Frequency")
                )
            })
        })
    })

    ## QC - DOWNLOAD BUTTON - Histogram (pdf)
    output$dlqchistpdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqchistpdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - DOWNLOAD PLOT - Histogram (pdf)
    output$dlqchistpdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-histogram.pdf")
        },
        content = function(file) {
            pdf(file, width = 9, height = 6.5, onefile = FALSE) # open the pdf device
            qcHist(
                tmp = ddstran()[[1]],
                lab = ddstran()[[2]],
                tran = input$transform
            )
            dev.off()
        }
    )

    ## QC - DOWNLOAD BUTTON - Histogram (png)
    output$dlqchistpng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqchistpngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - DOWNLOAD PLOT - Histogram (png)
    output$dlqchistpngimg <- downloadHandler(
        filename =  function() {
            paste("qc-histogram.png")
        },
        content = function(file) {
            png(file, width = 900, height = 750)
            qcHist(
                tmp = ddstran()[[1]],
                lab = ddstran()[[2]],
                tran = input$transform
            )
            dev.off()
        }
    )

    ## QC - header (2) - barplot
    output$counttotal <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Total reads")
        }
    })

    ## QC - visualize - barplot
    output$barplot <- renderPlotly({
        withProgress(message = "Compiling barplots...", value = 0, {
            incProgress()
            input$goqc
            tmp <- ddstran()[[1]]
            tmp <- assay(tmp)
            lab <- ddstran()[[2]]
            isolate({
                incProgress()
                dds <- ddsout()[[1]]
                bar <- as.data.frame(assay(dds))
                bar <- colSums(bar)
                bar <- as.data.frame(t(bar))
                bar <- gather(bar)
                plot_ly(
                    bar,
                    type = "bar",
                    y = ~value,
                    x = ~key,
                    color = ~key
                ) %>%
                layout(
                    margin = list(b = 90),
                    xaxis = list(title = "", tickangle = -45),
                    yaxis = list(title = "Counts")
                )
            })
        })
    })

    ## QC - DOWNLOAD BUTTON - Barplot (pdf)
    output$dlqcbarplotpdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcbarplotpdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - DOWNLOAD PLOT - Barplot (pdf)
    output$dlqcbarplotpdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-barplot.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6.5, onefile = FALSE) # open the pdf device
            qcBarplot(
                tmp = ddsout()[[1]]
            )
            dev.off()
        }
    )

    ## QC - DOWNLOAD BUTTON - Barplot (png)
    output$dlqcbarplotpng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcbarplotpngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - DOWNLOAD PLOT - Barplot (png)
    output$dlqcbarplotpngimg <- downloadHandler(
        filename =  function() {
            paste("qc-barplot.png")
        },
        content = function(file) {
            png(file, width = 900, height = 750)
            qcBarplot(
                tmp = ddsout()[[1]]
            )
            dev.off()
        }
    )

    ## QC - header (2) - barplot
    output$counttotal <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Total reads")
        }
    })

    ## QC - header (2) - PCA
    output$headpca <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Principal Component Analysis")
        }
    })

    ## QC - select input - choose factor - PCA
    output$pcafact <- renderUI({
        tmp <- ddsout()[[2]]
        selectInput(
            inputId = "pcafact",
            label = "Choose factor",
            choices = colnames(tmp)
        )
    })

    ## QC - visualize - PCA
    output$pca <- renderPlotly({
        validate(
            need(input$pcafact != "", "")
        )
        tmp <- ddstran()[[1]]
        lab <- ddstran()[[2]]
        validate(
            need(
                expr = class(tmp) == "DESeqTransform",
                message = "Please transform raw counts to visualize PCA."
            )
        )
        pca.lab <- plotPCA(
            tmp,
            intgroup = input$pcafact
        )
        pca <- plotPCA(
            tmp,
            intgroup = input$pcafact,
            returnData = TRUE
        )
        tooltips <- paste0(
            "<b>Sample:</b> ", rownames(pca), "<br />",
            "<b>PC1:</b> ", round(pca$PC1, 3), "<br />",
            "<b>PC2:</b> ", round(pca$PC2, 3)
        )
        plot_ly(
            data = pca,
            type = "scatter",
            mode = "markers",
            x = ~PC1,
            y = ~PC2,
            symbol = ~group,
            marker = list(size = 9),
            text = tooltips,
            hoverinfo = "text"
        ) %>%
        layout(
            xaxis = list(title = pca.lab$labels$x),
            yaxis = list(title = pca.lab$labels$y)
        )
    })

    ## QC - Show download button - PCA (PDF)
    output$dlqcpcapdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcpcapdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - Download plot - PCA (PDF)
    output$dlqcpcapdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-pca.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
            qcPCAPlot(
                tmp = ddstran()[[1]],
                pcafact = input$pcafact
            )
            dev.off()
        }
    )

    ## QC - Show download button - PCA (PNG)
    output$dlqcpcapng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcpcapngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - Download plot - PCA (PNG)
    output$dlqcpcapngimg <- downloadHandler(
        filename =  function() {
            paste("qc-pca.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            qcPCAPlot(
                tmp = ddstran()[[1]],
                pcafact = input$pcafact
            )
            dev.off()
        }
    )

    ## QC - header (2) - MDS
    output$headmds <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Multidimensional Scaling")
        }
    })

    ## QC - select input - choose factor - MDS
    output$mdsfact <- renderUI({
        tmp <- ddsout()[[2]]
        selectInput(
            inputId = "mdsfact",
            label = "Choose factor",
            choices = colnames(tmp)
        )
    })

    ## QC - visualize - MDS
    output$mds <- renderPlotly({
        tmp <- ddstran()[[1]]
        lab <- ddstran()[[2]]
        validate(
            need(
                expr = class(tmp) == "DESeqTransform",
                message = "Please transform raw counts to visualize MDS."
            )
        )
        mds <- dist(t(assay(tmp)))
        mds <- as.matrix(mds)
        mds <- as.data.frame(colData(tmp)) %>%
            cbind(cmdscale(mds))

        tooltips <- paste0(
            "<b>Sample:</b> ", rownames(mds), "<br />",
            "<b>Coord. 1:</b> ", round(mds[, "1"], 3), "<br />",
            "<b>Coord. 2:</b> ", round(mds[, "2"], 3)
        )
        plot_ly(
            data = mds,
            type = "scatter",
            mode = "markers",
            x = mds[, "1"],
            y = mds[, "2"],
            symbol = mds[, input$mdsfact],
            marker = list(size = 9),
            text = tooltips,
            hoverinfo = "text"
        ) %>%
        layout(
            xaxis = list(title = "MDS coordinate 1"),
            yaxis = list(title = "MDS coordinate 2")
        )
    })

    ## QC - Show download button - MDS (PDF)
    output$dlqcmdspdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcmdspdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - Download plot - MDS (PDF)
    output$dlqcmdspdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-mds.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
            qcMDSPlot(
                tmp = ddstran()[[1]],
                mdsfact = input$mdsfact
            )
            dev.off()
        }
    )



    ## QC - Show download button - MDS (PNG)
    output$dlqcmdspng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcmdspngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - Download plot - MDS (PNG)
    output$dlqcmdspngimg <- downloadHandler(
        filename =  function() {
            paste("qc-mds.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            qcMDSPlot(
                tmp = ddstran()[[1]],
                mdsfact = input$mdsfact
            )
            dev.off()
        }
    )


    ## QC - header (2) - tSNE
    output$headTSNE <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("T-Distributed Stochastic Neighbor Embedding")
        }
    })

    ## QC - select input - choose factor - tSNE
    output$tsnefact <- renderUI({
      tmp <- ddsout()[[2]]
      selectInput(
        inputId = "tsnefact",
        label = "Choose factor",
        choices = colnames(tmp)
      )
    })

    # QC - select dimensionality - tSNE
    output$tsneDim <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            radioButtons(
                inputId = "tsneDim",
                label = "Choose dimensions",
                choices = c(2, 3),
                inline = TRUE
            )
        }
    })

    ## QC - select perplexity - tSNE
    output$tsnePerp <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            textInput(
                inputId = "tsnePerp",
                label = "Enter perplexity",
                value = 1,
                width = "100%"
            )
        }
    })

    # QC - visualize - tSNE
    tsneout <- reactive({
        if (input$goqc == 0) {
            return()
        } else {
            tmp <- ddstran()[[1]]
            lab <- ddstran()[[2]]
            coldata <- colData(tmp)
            validate(
                need(
                    expr = class(tmp) == "DESeqTransform",
                    message = "Please transform raw counts to visualize tSNE."
                )
            )
            tsne <- as.matrix(assay(tmp))
            tsne <- dist(t(tsne))
            req(input$tsnePerp)
            tmp1 <- nrow(coldata) - 1
            tmp2 <- 3 * as.numeric(input$tsnePerp)
            if (tmp1 < tmp2) {
                tsneOut <- "na"
            } else {
                ### Maybe change verbose to TRUE later...
                tsneOut <- Rtsne::Rtsne(
                    X = tsne,
                    dims = 3,
                    perplexity = as.numeric(input$tsnePerp),
                    is_distance = TRUE,
                    verbose = FALSE,
                    pca = FALSE
                )
                tsneOut <- tsneOut$Y
            }

            return(list(tsneOut))
        }
    })

    output$tsnePerpCheck <- renderUI({
        tsneOut <- tsneout()[[1]]

        if (tsneOut == "na") {
            p(
                em(
                    paste0(
                        "Your perplexity value is too large. Please ",
                        "try another value that is lower!"
                    )
                ),
                style = "color:grey"
            )
        } else {
            return()
        }
    })

    output$tsnePlot <- renderPlotly({
        req(tsneout())
        tsneOut <- tsneout()[[1]]

        if (tsneOut == "na") {
            return()
        } else {
            req(input$tsneDim)
            if (as.numeric(input$tsneDim) == 2) {
                tmp <- ddstran()[[1]]
                lab <- ddstran()[[2]]
                coldata <- colData(tmp)

                colors = grDevices::rainbow(
                    length(unique(coldata[, input$tsnefact]))
                )
                names(colors) = unique(coldata[, input$tsnefact])

                tsneOutDF <- as.data.frame(tsneOut)
                plot_ly(
                    data = tsneOutDF,
                    type = "scatter",
                    mode = "markers",
                    x = tsneOutDF$V1,
                    y = tsneOutDF$V2,
                    symbol = names(colors[coldata[, input$tsnefact]]),
                    marker = list(size = 9)
                    # text = tooltips,
                    # hoverinfo = "text"
                ) %>%
                    layout(
                        xaxis = list(title = "tSNE coordinate 1"),
                        yaxis = list(title = "tSNE coordinate 2")
                    )
            } else if (as.numeric(input$tsneDim) == 3) {
                tmp <- ddstran()[[1]]
                lab <- ddstran()[[2]]
                coldata <- colData(tmp)

                colors = grDevices::rainbow(
                    length(unique(coldata[, input$tsnefact]))
                )
                names(colors) = unique(coldata[, input$tsnefact])

                tsneOutDF <- as.data.frame(tsneOut)
                plot_ly(
                    data = tsneOutDF,
                    mode = "markers",
                    x = tsneOutDF$V1,
                    y = tsneOutDF$V2,
                    z = tsneOutDF$V3,
                    symbol = names(colors[coldata[, input$tsnefact]]),
                    marker = list(size = 9)
                    # text = tooltips,
                    # hoverinfo = "text"
                ) %>%
                    add_markers() %>%
                    layout(
                        scene = list(
                            xaxis = list(title = "tSNE coordinate 1"),
                            yaxis = list(title = "tSNE coordinate 2"),
                            zaxis = list(title = "tSNE coordinate 3")
                        )
                    )
            }
        }
    })


    ## QC - Show download button - tSNE (PDF)
    output$dlqctsnepdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqctsnepdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## QC - Download plot - tSNE (PDF)
    output$dlqctsnepdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-tsne.pdf")
        },
        content = function(file) {
            tsneOut <- tsneout()[[1]]
            tmp <- ddstran()[[1]]
            lab <- ddstran()[[2]]
            coldata <- colData(tmp)

            colors = grDevices::rainbow(
                length(unique(coldata[, input$tsnefact]))
            )
            names(colors) = unique(coldata[, input$tsnefact])
            tsneOutDF <- as.data.frame(tsneOut)

            pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
            qcTSNEPlot(
                tsneDF = tsneOutDF,
                col = colors,
                coldata = coldata,
                tsnefact = input$tsnefact
            )
            dev.off()
        }
    )



    ## QC - Show download button - tsne (PNG)
    output$dlqctsnepng <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqctsnepngimg", "Download Static R Plot (PNG)")
        }
    })

    ## QC - Download plot - tsne (PNG)
    output$dlqctsnepngimg <- downloadHandler(
        filename =  function() {
            paste("qc-tsne.png")
        },
        content = function(file) {
            tsneOut <- tsneout()[[1]]
            tmp <- ddstran()[[1]]
            lab <- ddstran()[[2]]
            coldata <- colData(tmp)

            colors = grDevices::rainbow(
                length(unique(coldata[, input$tsnefact]))
            )
            names(colors) = unique(coldata[, input$tsnefact])
            tsneOutDF <- as.data.frame(tsneOut)

            png(file, width = 700, height = 650)
            p <- qcTSNEPlot(
                tsneDF = tsneOutDF,
                col = colors,
                coldata = coldata,
                tsnefact = input$tsnefact
            )
            print(p)
            dev.off()
        }
    )

    output$tsneDebug <- renderPrint({
        tsneout()[[1]]
    })





    ###################################################################
    ###################################################################
    ### SECTION 03 - DIFFERENTIAL GENE EXPRESSION ANALYSIS (DGE)
    ###################################################################
    ###################################################################

    ## DGE - Choose analytical methods
    output$dgemethod <- renderUI({
        validate(
            need(input$dgeexpsetup != "", "")
        )
        if (input$dgeexpsetup != "exp5" & input$dgeexpsetup != "exp6") {
            selectInput(
                inputId = "dgemethod",
                label = "Choose method",
                choices = c(
                    "DESeq2" = "deseq",
                    "edgeR" = "edger",
                    "limma-voom" = "limma"
                ),
                selected = "deseq"
            )
        } else {
            selectInput(
                inputId = "dgemethod",
                label = "Choose method",
                choices = c(
                    "DESeq2" = "deseq",
                    "edgeR" = "edger"
                ),
                selected = "deseq"
            )
        }
    })

    ## DGE - exp. setup 1 - two group comparisons - factor choice
    output$dgeexp1a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp1") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp1a",
                    label = "Choose factor",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 1 - two group comparisons - choose comparisons
    output$dgeexp1b <- renderUI({
        validate(
            need(input$dgeexp1a != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp1") {
                perm <- ddsout()[[2]]
                perm <- as.vector(unique(perm[, input$dgeexp1a]))
                perm <- permutations(n = length(perm), r = 2, v = perm)
                perm <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")
                checkboxGroupInput(
                    inputId = "dgeexp1b",
                    label = "Choose comparisons you want made",
                    choices = as.list(perm)
                )
            }
        }
    })

    ## DGE - exp. setup 2 - mult. group comparisons - factor choice A
    output$dgeexp2a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
                p(
                    em(
                        "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
                    ),
                    style = "color:grey"
                )
            } else if (input$dgeexpsetup == "exp2") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp2a",
                    label = "Choose factor (A)",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 2 - mult. group comparisons - factor choice B
    output$dgeexp2b <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp2") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp2b",
                    label = "Choose factor (B)",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 2 - mult. group comparisons - group comb. choices
    output$dgeexp2c <- renderUI({
        validate(
            need(input$dgeexp2a != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp2") {
                perm <- ddsout()[[2]]
                group.c <- factor(
                    paste(perm[, input$dgeexp2a], perm[, input$dgeexp2b], sep = "_")
                )
                perm <- levels(group.c)
                perm <- permutations(n = length(perm), r = 2, v = perm)
                perm <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")
                checkboxGroupInput(
                    inputId = "dgeexp2c",
                    label = "Choose comparisons you want made",
                    choices = as.list(perm)
                )
            }
        }
    })

    ## DGE - exp. setup 3 - interaction - factor choice A
    output$dgeexp3a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
                p(
                    em(
                        "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
                    ),
                    style = "color:grey"
                )
            } else if (input$dgeexpsetup == "exp3") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp3a",
                    label = "Choose factor (A)",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 3 - interaction - factor choice B
    output$dgeexp3b <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp3") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp3b",
                    label = "Choose factor (B)",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 3 - interaction - reference level for factor A
    output$dgeexp3c <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp3") {
                tmp <- ddsout()[[2]]
                tmp <- unique(tmp[, input$dgeexp3a])
                selectInput(
                    inputId = "dgeexp3c",
                    label = "Choose reference level for factor A",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 3 - interaction - reference level for factor B
    output$dgeexp3d <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp3") {
                tmp <- ddsout()[[2]]
                tmp <- unique(tmp[, input$dgeexp3b])
                selectInput(
                    inputId = "dgeexp3d",
                    label = "Choose reference level for factor B",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 4 - additive model - blocking factor
    output$dgeexp4a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
                p(
                    em(
                        "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
                    ),
                    style = "color:grey"
                )
            } else if (input$dgeexpsetup == "exp4") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp4a",
                    label = "Choose your blocking or subjet factor",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 4 - additive model - treatment factor
    output$dgeexp4b <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp4") {
                tmp <- ddsout()[[2]]
                selectInput(
                    inputId = "dgeexp4b",
                    label = "Choose your treatment factor",
                    choices = colnames(tmp)
                )
            }
        }
    })

    ## DGE - exp. setup 4 - additive model - reference level for blocking factor
    output$dgeexp4c <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp4") {
                tmp <- ddsout()[[2]]
                tmp <- unique(tmp[, input$dgeexp4a])
                selectInput(
                    inputId = "dgeexp4c",
                    label = "Choose reference level for blocking factor",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 4 - additive model - reference level for treatment factor
    output$dgeexp4d <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp4") {
                tmp <- ddsout()[[2]]
                tmp <- unique(tmp[, input$dgeexp4b])
                selectInput(
                    inputId = "dgeexp4d",
                    label = "Choose reference level for treatment factor",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 5 - main effect (ME) - choose ME
    output$dgeexp5a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp5") {
                tmp <- ddsout()[[2]]
                tmp <- unique(colnames(tmp))
                selectInput(
                    inputId = "dgeexp5a",
                    label = "Choose main effect",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 5 - ME - choose ME reference
    output$dgeexp5b <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp5") {
                tmp <- ddsout()[[2]]
                tmp <- levels(tmp[, input$dgeexp5a])
                selectInput(
                    inputId = "dgeexp5b",
                    label = "Choose main effect reference level",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 5 - ME - choose contrasts
    output$dgeexp5c <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp5") {
                tmp1 <- ddsout()[[2]]
                tmp2 <- input$dgeexp5a
                tmp3 <- levels(tmp1[, tmp2])
                tmp3 <- tmp3[which(tmp3 != input$dgeexp5b)]
                tmp3 <- paste0(tmp3, "_VS_", input$dgeexp5b)
                checkboxGroupInput(
                    inputId = "dgeexp5c",
                    label = "Choose comparisons you want made",
                    choices = as.list(tmp3)
                )
            }
        }
    })

    ## DGE - exp. setup 6 - ME + group fact - choose ME
    output$dgeexp6a <- renderUI({
        validate(
            need(input$dgemethod != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp6") {
                tmp <- ddsout()[[2]]
                tmp <- unique(colnames(tmp))
                selectInput(
                    inputId = "dgeexp6a",
                    label = "Choose grouping factor",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 6 - ME + group fact - choose group fact
    output$dgeexp6b <- renderUI({
        validate(
            need(input$dgemethod != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp6") {
                tmp <- ddsout()[[2]]
                tmp <- levels(tmp[, input$dgeexp6a])
                selectInput(
                    inputId = "dgeexp6b",
                    label = "Choose grouping factor level",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 6 - ME + group fact - choose ME reference
    output$dgeexp6c <- renderUI({
        validate(
            need(input$dgemethod != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                tmp <- ddsout()[[2]]
                tmp <- unique(colnames(tmp))
                selectInput(
                    inputId = "dgeexp6c",
                    label = "Choose main effect",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 6 - ME + group fact - choose group fact level
    output$dgeexp6d <- renderUI({
        validate(
            need(
                expr = !is.null(input$dgeexp6a),
                message = ""
            )
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp6") {
                tmp <- ddsout()[[2]]
                tmp <- tmp[which(tmp[, input$dgeexp6a] == input$dgeexp6b), ]
                tmp[] <- lapply(tmp, function(x) if(is.factor(x)) factor(x) else x)
                tmp <- levels(tmp[, input$dgeexp6c])
                selectInput(
                    inputId = "dgeexp6d",
                    label = "Choose main effect reference level",
                    choices = tmp
                )
            }
        }
    })

    ## DGE - exp. setup 6 - ME group fact - choose contrasts
    output$dgeexp6e <- renderUI({
        validate(
            need(
                expr = !is.null(input$dgeexp6a),
                message = ""
            )
        )
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                tmp1 <- ddsout()[[2]]
                tmp2 <- tmp1[which(tmp1[, input$dgeexp6a] == input$dgeexp6b), ]
                tmp2[] <- lapply(tmp2, function(x) if(is.factor(x)) factor(x) else x)
                tmp2 <- levels(tmp2[, input$dgeexp6c])
                tmp2 <- tmp2[which(tmp2 != input$dgeexp6d)]
                tmp3 <- paste0(tmp2, "_VS_", input$dgeexp6d)
                checkboxGroupInput(
                    inputId = "dgeexp6e",
                    label = "Choose comparisons you want made",
                    choices = as.list(tmp3)
                )
            }
        }
    })

    ## DGE - exp. setup 7 - user input
    output$dgeexp7a <- renderUI({
        validate(
            need(input$dgemethod != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else if (input$dgeexpsetup != "exp7") {
            return()
        } else {
            fileInput(
                inputId = "mod.matrix",
                label = "Submit model matrix (CSV)",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            )
        }
    })

    ## DGE - exp. setup - formula - header
    output$dgeexpformhead <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup != "exp1" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp1" & ncol(ddsout()[[2]] < 2)){
                h5(
                    strong("Your linear model will look like this:")
                )
            } else if (input$dgeexpsetup == "exp7") {
                return()
            } else {
                h5(
                    strong("Your linear model will look like this:")
                )
            }
        }
    })

    ## DGE - exp. setup - formula - formula 1
    output$dgeexpform1 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp1") {
                code(
                    paste0(" ~ ", input$dgeexp1a)
                )
            }
        }
    })

    ## DGE - exp. setup - formula - formula 2
    output$dgeexpform2 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp2") {
                code(
                    paste0(" ~ ", input$dgeexp2a, "_", input$dgeexp2b)
                )
            }
        }
    })

    ## DGE - exp. setup 3 - interaction - formula layout
    output$dgeexpform3 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp3") {
                code(
                    paste0(
                        " ~ ", input$dgeexp3a,
                        " + ", input$dgeexp3b,
                        " + ", input$dgeexp3a,
                        ":", input$dgeexp3b
                    )
                )
            }
        }
    })

    ## DGE - exp. setup 4 - added effects - formula layout
    output$dgeexpform4 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
                return()
            } else if (input$dgeexpsetup == "exp4") {
                code(
                    paste0(" ~ ", input$dgeexp4a, " + ", input$dgeexp4b)
                )
            }
        }
    })

    ## DGE - exp. setup 5 - ME - formula layout
    output$dgeexpform5 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp5") {
                code(
                    paste0(" ~ ", input$dgeexp5a, " (as main effect)")
                )
            }
        }
    })

    ## DGE - exp. setup 6A - ME + group fact. - formula layout
    output$dgeexpform6a <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                code(
                    paste0(
                        " ~ ", input$dgeexp6c, " (as main effect)"
                    )
                )
            }
        }
    })

    ## DGE - exp. setup 6B - ME + group fact. - formula layout
    output$dgeexpform6b <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                code(
                    paste0(
                        "* Limited by: "
                    )
                )
            }
        }
    })

    ## DGE - exp. setup 6C - ME + group fact. - formula layout
    output$dgeexpform6c <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                code(
                    paste(
                        "...factor: ",
                        input$dgeexp6a
                    )
                )
            }
        }
    })

    ## DGE - exp. setup 6D - ME + group fact. - formula layout
    output$dgeexpform6d <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            if (input$dgeexpsetup == "exp6") {
                code(
                    paste(
                        "...factor level: ", input$dgeexp6b
                    )
                )
            }
        }
    })


    ## DGE - edgeR normalization option
    output$dgeexpedgernorm <- renderUI({
        validate(
            need(input$dgemethod != "", "")
        )
        if (input$dgemethod != "edger") {
            return()
        } else {
            selectInput(
                inputId = "dgeexpedgernorm",
                label = "Choose normalization type",
                choices = c(
                    "TMM" = "TMM",
                    "RLE" = "RLE",
                    "upperquartile" = "upperquartile",
                    "none" = "none"
                )
            )
        }
    })

    ## DGE - Choose contrasts
    output$dgemaincontrasts <- renderUI({
        if (input$godge == 0) {
            return()
        } else {
            tmp <- dgeout1()[[2]]
            tmp <- colnames(tmp)
            selectInput(
                inputId = "dgemaincontrasts",
                label = "Choose contrast",
                choices = tmp,
                width = "400px"
            )
        }
    })

    ## DGE - user input for model matrix
    mod.matrix <- eventReactive(input$godge, {
        mod.matrix <- input$mod.matrix
        mod.matrix <- as.matrix(
            read.csv(
                mod.matrix$datapath, header = TRUE, row.names = 1
            )
        )
        return(list(mod.matrix))
    })

    # ## DEBUG ##
    # output$debugdge2 <- renderPrint({
    #   if (input$godge == 0) {
    #     return()
    #   } else {
    #     mod.matrix()[[1]]
    #   }
    # })

    ## DGE - analysis - reactive
    dgeout1 <- eventReactive(input$godge, {
        cts <- ddsout()[[3]]
        coldata <- ddsout()[[2]]
        if (input$dgemethod == "limma") {
            if (input$dgeexpsetup == "exp1") {
                withProgress(message = "Running limma-voom...", value = 0, {
                    incProgress(1/2)
                    de.genes <- limma.exp1(
                        fact = input$dgeexp1a,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp1b
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp2") {
                withProgress(message = "Running limma-voom...", value = 0, {
                    incProgress(1/2)
                    de.genes <- limma.exp2(
                        fact1 = input$dgeexp2a,
                        fact2 = input$dgeexp2b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp2c
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp3") {
                withProgress(message = "Running limma-voom...", value = 0, {
                    incProgress(1/2)
                    de.genes <- limma.exp3(
                        fact1 = input$dgeexp3a,
                        fact2 = input$dgeexp3b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp3c,
                        fact2.rlvl = input$dgeexp3d
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp4") {
                withProgress(message = "Running limma-voom...", value = 0, {
                    incProgress(1/2)
                    de.genes <- limma.exp4(
                        fact1 = input$dgeexp4a,
                        fact2 = input$dgeexp4b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp4c,
                        fact2.rlvl = input$dgeexp4d
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp7") {
                withProgress(message = "Running limma-voom...", value = 0, {
                    incProgress(1/2)
                    de.genes <- limma.exp7(
                        cts = cts,
                        mod.matrix = mod.matrix()[[1]]
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            }
        } else if (input$dgemethod == "edger") {
            if (input$dgeexpsetup == "exp1") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp1(
                        fact = input$dgeexp1a,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp1b,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp2") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp2(
                        fact1 = input$dgeexp2a,
                        fact2 = input$dgeexp2b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp2c,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp3") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp3(
                        fact1 = input$dgeexp3a,
                        fact2 = input$dgeexp3b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp3c,
                        fact2.rlvl = input$dgeexp3d,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp4") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp4(
                        fact1 = input$dgeexp4a,
                        fact2 = input$dgeexp4b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp4c,
                        fact2.rlvl = input$dgeexp4d,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp5") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp5(
                        fact = input$dgeexp5a,
                        fact.levl = input$dgeexp5b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp5c,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp6") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp6(
                        me.fact = input$dgeexp6c,
                        me.levl = input$dgeexp6d,
                        gp.fact = input$dgeexp6a,
                        gp.levl = input$dgeexp6b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp6e,
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp7") {
                withProgress(message = "Running edgeR...", value = 0, {
                    incProgress(1/2)
                    de.genes <- edger.exp7(
                        cts = cts,
                        mod.matrix = mod.matrix()[[1]],
                        norm = input$dgeexpedgernorm
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            }
        } else if (input$dgemethod == "deseq") {
            if (input$dgeexpsetup == "exp1") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp1(
                        fact = input$dgeexp1a,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp1b
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp2") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp2(
                        fact1 = input$dgeexp2a,
                        fact2 = input$dgeexp2b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp2c
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp3") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp3(
                        fact1 = input$dgeexp3a,
                        fact2 = input$dgeexp3b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp3c,
                        fact2.rlvl = input$dgeexp3d
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp4") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp4(
                        fact1 = input$dgeexp4a,
                        fact2 = input$dgeexp4b,
                        cts = cts,
                        coldata = coldata,
                        fact1.rlvl = input$dgeexp4c,
                        fact2.rlvl = input$dgeexp4d
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp5") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp5(
                        fact = input$dgeexp5a,
                        fact.levl = input$dgeexp5b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp5c
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp6") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp6(
                        me.fact = input$dgeexp6c,
                        me.levl = input$dgeexp6d,
                        gp.fact = input$dgeexp6a,
                        gp.levl = input$dgeexp6b,
                        cts = cts,
                        coldata = coldata,
                        perm.h = input$dgeexp6e
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            } else if (input$dgeexpsetup == "exp7") {
                withProgress(message = "Running DESeq2...", value = 0, {
                    incProgress(1/2)
                    de.genes <- deseq.exp7(
                        cts = cts,
                        coldata = coldata,
                        mod.matrix = mod.matrix()[[1]]
                    )
                    fit.names <- de.genes[[2]]
                    fit.cont <- de.genes[[1]]
                    incProgress(2/2)
                })
            }
        }
        return(list(fit.cont, fit.names))
    })

    ## DGE - header (2) - DGE Overview
    output$headdgeoverview <- renderUI({
        if(input$godge == 0) {
            p(
                br(),
                em(
                    "Choose an experimental setup with parameters and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("DGE Overview")
        }
    })

    ## DGE - header (2) - DGE Overview
    output$headdgeplots <- renderUI({
        if(input$godge == 0) {
            p(
                br(),
                em(
                    "Choose an experimental setup with parameters and click the 'submit' button to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Interactive Plots")
        }
    })

    ## DGE - contrast table
    dgeout2 <- reactive({
        if(input$godge == 0) {
            return()
        } else {
            isolate({
                expset <- input$dgeexpsetup
            })
            contTable <- getContTable(
                de.genes = dgeout1()[[1]],
                coef = input$dgemaincontrasts,
                cts = ddsout()[[3]],
                expset = expset,
                design = dgeout1()[[2]],
                fact = input$dgeexp1a,
                fact5 = input$dgeexp5a,
                fact6 = input$dgeexp6c
            )
            return(list(contTable))
        }
    })

    ## DGE - select input - Visualization type
    output$vistype <- renderUI({
        if(input$godge == 0) {
            return()
        } else {
            radioButtons(
                inputId = "plottype",
                label = "Choose plot type",
                choices = c(
                    "MA plot" = "maplot",
                    "Volcano plot" = "volplot"
                ),
                selected = "maplot",
                inline = TRUE
            )
        }
    })

    ## DGE - Filter Data
    dgeout3 <- reactive({
        if (input$godge == 0) {
            return()
        } else {
            tmp <- dgeout2()[[1]]
            # is.na(tmp$padj) <- 1
            # tmp$padj <- round(tmp$padj, 3)
            padj <- as.numeric(input$dgepadjcutoff)
            lfc <- as.numeric(input$dgefcmin)
            tmp <- tmp[abs(tmp$log2FoldChange) >= lfc, ]
            tmp <- tmp[tmp$padj <= padj, ]
            return(tmp)
        }
    })

    ## DGE - Shared data
    share <- reactive({
        SharedData$new(dgeout3())
    })

    ## DGE - Conditional plot
    output$dgeplot <- renderPlotly({
        validate(
            need(input$dgemaincontrasts != "", "")
        )
        check <- dgeout3()
        check <- nrow(check)
        validate(
            need(
                expr = check > 1,
                message = paste(
                    "Note:",
                    "It seems that you have no differentially",
                    "expressed IDs. There are zero differentially",
                    "expressed points to plot. You may still download",
                    "static images of the entire dataset as either",
                    "PDF or PNG files."
                )
            )
        )
        if (input$godge == 0) {
            return()
        } else {
            s <- input$mytable_rows_selected
            datobj <- dgeout3()
            tooltips <- paste0(
                "<b>ID:</b> ", datobj$id, "<br />",
                "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
                "<b>PADJ:</b> ", round(datobj$padj, 3), "<br />",
                "<b>BM:</b> ", round(datobj$baseMean, 3)
            )
            if (input$godge == 0) {
                return()
            } else {
                if (input$plottype == "maplot") {
                    if (!length(s)) {
                        p <- share() %>%
                            plot_ly(
                                x = ~log10(baseMean),
                                y = ~log2FoldChange,
                                type = "scatter",
                                mode = "markers",
                                color = I("darkgray"),
                                text = tooltips,
                                hoverinfo = "text",
                                name = "Unfiltered") %>%
                            layout(
                                showlegend = TRUE,
                                xaxis = list(title = "log<sub>10</sub>(baseMean)"),
                                yaxis = list(title = "log<sub>2</sub>(fold change)")
                            ) %>%
                            highlight(
                                "plotly_selected",
                                color = I("royalblue1"),
                                selected = attrs_selected(
                                    name = "Filtered"
                                )
                            )
                    } else if (length(s)) {
                        pp <- dgeout3() %>%
                            plot_ly() %>%
                            add_trace(
                                x = ~log10(baseMean),
                                y = ~log2FoldChange,
                                type = "scatter",
                                mode = "markers",
                                color = I("darkgray"),
                                text = tooltips,
                                hoverinfo = "text",
                                name = "Unfiltered") %>%
                            layout(
                                showlegend = TRUE,
                                xaxis = list(title = "log<sub>10</sub>(baseMean)"),
                                yaxis = list(title = "log<sub>2</sub>(fold change)")
                            )

                        # selected data
                        pp <- add_trace(
                            pp,
                            data = dgeout3()[s, , drop = FALSE],
                            x = ~log10(baseMean),
                            y = ~log2FoldChange,
                            type = "scatter",
                            mode = "markers",
                            color = I("royalblue1"),
                            # text = tooltips,
                            size = I(8),
                            name = "Filtered"
                        )
                    }
                } else if (input$plottype == "volplot") {
                    if (!length(s)) {
                        p <- share() %>%
                            plot_ly(
                                x = ~log2FoldChange,
                                y = ~-log10(pvalue),
                                type = "scatter",
                                mode = "markers",
                                color = I("darkgray"),
                                text = tooltips,
                                hoverinfo = "text",
                                name = "Unfiltered") %>%
                            layout(
                                showlegend = TRUE,
                                xaxis = list(title = "log<sub>2</sub>(fold change)"),
                                yaxis = list(title = "-log<sub>10</sub>(p-value)")
                            ) %>%
                            highlight(
                                "plotly_selected",
                                color = I("royalblue1"),
                                selected = attrs_selected(
                                    name = "Filtered"
                                )
                            )
                    } else if (length(s)) {
                        pp <- dgeout3() %>%
                            plot_ly() %>%
                            add_trace(
                                x = ~log2FoldChange,
                                y = ~-log10(pvalue),
                                type = "scatter",
                                mode = "markers",
                                color = I("darkgray"),
                                text = tooltips,
                                hoverinfo = "text",
                                name = "Unfiltered") %>%
                            layout(
                                showlegend = TRUE,
                                xaxis = list(title = "log<sub>2</sub>(fold change)"),
                                yaxis = list(title = "-log<sub>10</sub>(p-value)")
                            )

                        # selected data
                        pp <- add_trace(
                            pp,
                            data = dgeout3()[s, , drop = FALSE],
                            x = ~log2FoldChange,
                            y = ~-log10(pvalue),
                            mode = "markers",
                            color = I("royalblue1"),
                            # text = tooltips,
                            size = I(8),
                            name = "Filtered"
                        )
                    }
                }
            }
        }
    })

    ## DGE - Conditional table
    output$mytable <- DT::renderDataTable({
        validate(
            need(input$dgemaincontrasts != "", "")
        )
        if (input$godge == 0) {
            return()
        } else {
            tmp <- dgeout3() %>% mutate_if(is.numeric, round, digits = 4)
            # tmp <- dgeout3()
            tmp <- tmp[complete.cases(tmp), ]
            m2 <- tmp[share()$selection(),]
            dt <- DT::datatable(tmp)
            if (NROW(m2) == 0) {
                dt
            } else {
                DT::formatStyle(
                    dt,
                    "id",
                    target = "row",
                    color = DT::styleEqual(m2$id, rep("white", length(m2$id))),
                    backgroundColor = DT::styleEqual(
                        m2$id,
                        rep("darkgray", length(m2$id))
                    )
                )
            }
        }
    })

    ## DGE - Show download button (FILTERED)
    output$downloadfilt <- renderUI({
        validate(
            need(
                expr = !is.null(input$godge),
                message = ""
            )
        )
        if(input$godge == 0) {
            return()
        } else {
            downloadButton("downfiltData", "Download Filtered Data")
        }
    })

    ## DGE - Download FILTERED data
    output$downfiltData <- downloadHandler(
        filename = function() {
            paste0(
                input$dgemaincontrasts,
                "filtered_results_",
                nrow(dgeout3()[input[["mytable_rows_all"]], ]),
                ".csv"
            )
        },
        content = function(file) {
            filtData <- dgeout3() %>% mutate_if(is.numeric, round, digits = 4)
            filtData <- filtData[complete.cases(filtData), ]
            write.csv(
                filtData,
                file,
                row.names = FALSE
            )
        }
    )

    ## DGE - Show download button (ALL)
    output$downloadall <- renderUI({
        validate(
            need(
                expr = !is.null(input$godge),
                message = ""
            )
        )
        if(input$godge == 0) {
            return()
        } else {
            downloadButton("downallData", "Download All Data")
        }
    })

    ## DGE - Download ALL data
    output$downallData <- downloadHandler(
        filename = function() {
            paste0(input$dgemaincontrasts, "_all_results.csv")
        },
        content = function(file) {
            write.csv(
                dgeout2()[[1]],
                file, row.names = FALSE
            )
        }
    )

    ## DGE - Generate overview data
    dgeover0 <- reactive({
        if (input$godge == 0) {
            return()
        } else{
            expset <- input$dgeexpsetup
            perm <- dgeout1()[[2]]
            perm <- colnames(perm)
            cont.ls <- list()
            for (i in perm) {
                cont.ls[[i]] <- contTable <- getContTable(
                    de.genes = dgeout1()[[1]],
                    coef = i,
                    cts = ddsout()[[3]],
                    expset = expset,
                    design = dgeout1()[[2]],
                    fact = input$dgeexp1a,
                    fact5 = input$dgeexp5a,
                    fact6 = input$dgeexp6c
                )
            }
            return(list(cont.ls))
        }
    })

    ## DGE - table generation
    dgeover <- reactive({
        p <- as.numeric(input$dgepadjcutoff)
        lf <- as.numeric(input$dgefcmin)
        test <- dgeOverTbl(
            cont.ls = dgeover0()[[1]],
            lf = lf,
            p = p
        )
        return(list(test))
    })

    ## DGE - visualization - DGE overview
    output$dgeplot2 <- renderPlotly({
        if (input$godge == 0) {
            return()
        } else {
            comp <- dgeover()[[1]]
            plot_ly(
                comp,
                type = "bar",
                y = ~value,
                x = ~contrast,
                color = ~variable
            ) %>%
            layout(
                margin = list(b = 300),
                xaxis = list(title = "", tickangle = -90),
                yaxis = list(title = "Number of IDs")
            )
        }
    })

    ## DGE - Show download button - DGE Overview (PDF)
    output$dldgeoverpdf <- renderUI({
        if(input$godge == 0) {
            return()
        } else {
            downloadButton("dldgeoverpdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## DGE - Download plot - DGE Overview (PDF)
    output$dldgeoverpdfimg <- downloadHandler(
        filename =  function() {
            paste("dge-overview.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6.5, onefile = FALSE) # open the pdf device
            dgeOverPlot(
                comp = dgeover()[[1]]
            )
            dev.off()
        }
    )

    ## DGE - Show download button - DGE Overview (PNG)
    output$dldgeoverpng <- renderUI({
        if(input$godge == 0) {
            return()
        } else {
            downloadButton("dldgeoverpngimg", "Download Static R Plot (PNG)")
        }
    })

    ## DGE - Download plot - DGE Overview (PNG)
    output$dldgeoverpngimg <- downloadHandler(
        filename =  function() {
            paste("dge-overview.png")
        },
        content = function(file) {
            png(file, width = 900, height = 750)
            dgeOverPlot(
                comp = dgeover()[[1]]
            )
            dev.off()
        }
    )

    ## DGE - Overview table
    output$dgeoverview <- renderUI({
        if (input$godge == 0) {
            return()
        } else {
            # cont <- dgeover()[[1]]
            # colnames(cont) <- c("Comparison", "Regulation", "IDs")
            DT::dataTableOutput("dgeoverviewtable")
        }
    })

    ## DGE - Overview table
    output$dgeoverviewtable <- DT::renderDataTable(DT::datatable({
        cont <- dgeover()[[1]]
        colnames(cont) <- c("Comparison", "Regulation", "IDs")
        cont
    }))

    ## DGE - Show download button (DGE Overview Table)
    output$dldgeoverviewtbl <- renderUI({
        if(input$godge == 0) {
            return()
        } else {
            downloadButton("dldgeoverviewtbl2", "Download Table")
        }
    })

    ## DGE - Download DGE Overview Table
    output$dldgeoverviewtbl2 <- downloadHandler(
        filename = function() {
            paste0("dge-overview.csv")
        },
        content = function(file) {
            cont <- dgeover()[[1]]
            colnames(cont) <- c("Comparison", "Regulation", "IDs")
            write.csv(
                cont,
                file,
                row.names = FALSE
            )
        }
    )

    ## DGE - Show Download Button - Conditional - MA or Volcano (PDF)
    output$dldgemavolpdf <- renderUI({
        validate(
            need(input$plottype != "", "")
        )
        if (input$plottype == "maplot") {
            ### DGE - Show download button - MA plot (PDF)
            if(input$goqc == 0) {
                return()
            } else {
                downloadButton("dldgemapdfimg", "Download MA Plot (PDF)")
            }
        } else if (input$plottype == "volplot") {
            ### DGE - Show download button - Volcano plot (PDF)
            if(input$goqc == 0) {
                return()
            } else {
                downloadButton("dldgevolpdfimg", "Download Volcano Plot (PDF)")
            }
        }
    })

    ## DGE - Show Download Button - Conditional - MA or Volcano (PNG)
    output$dldgemavolpng <- renderUI({
        validate(
            need(input$plottype != "", "")
        )
        if (input$plottype == "maplot") {
            ### DGE - Show download button - MA plot (PNG)
            if(input$goqc == 0) {
                return()
            } else {
                downloadButton("dldgemapngimg", "Download MA Plot (PNG)")
            }
        } else if (input$plottype == "volplot") {
            ### DGE - Show download button - Volcano plot (PNG)
            if(input$goqc == 0) {
                return()
            } else {
                downloadButton("dldgevolpngimg", "Download Volcano Plot (PNG)")
            }
        }
    })

    ## DGE - Download plot - MA plot (PDF)
    output$dldgemapdfimg <- downloadHandler(
        filename =  function() {
            paste("dge-ma-plot.pdf")
        },
        content = function(file) {
            pdf(file, width = 9, height = 6.5, onefile = FALSE)
            dgeMAPlot(
                dgeout2 = dgeout2()[[1]],
                p = as.numeric(input$dgepadjcutoff),
                l = as.numeric(input$dgefcmin),
                cont = input$dgemaincontrasts
            )
            dev.off()
        }
    )

    ## DGE - Download plot - MA plot (PNG)
    output$dldgemapngimg <- downloadHandler(
        filename =  function() {
            paste("dge-ma-plot.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            dgeMAPlot(
                dgeout2 = dgeout2()[[1]],
                p = as.numeric(input$dgepadjcutoff),
                l = as.numeric(input$dgefcmin),
                cont = input$dgemaincontrasts
            )
            dev.off()
        }
    )

    ## DGE - Download plot - Volcano plot (PDF)
    output$dldgevolpdfimg <- downloadHandler(
        filename =  function() {
            paste("dge-vol-plot.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6.5, onefile = FALSE)
            dgeVolPlot(
                dgeout2 = dgeout2()[[1]],
                p = as.numeric(input$dgepadjcutoff),
                l = as.numeric(input$dgefcmin),
                cont = input$dgemaincontrasts
            )
            dev.off()
        }
    )

    ## DGE - Download plot - Volcano plot (PNG)
    output$dldgevolpngimg <- downloadHandler(
        filename =  function() {
            paste("dge-vol-plot.png")
        },
        content = function(file) {
            png(file, width = 800, height = 650)
            dgeVolPlot(
                dgeout2 = dgeout2()[[1]],
                p = as.numeric(input$dgepadjcutoff),
                l = as.numeric(input$dgefcmin),
                cont = input$dgemaincontrasts
            )
            dev.off()
        }
    )





    ###################################################################
    ###################################################################
    ### SECTION 04 - HEATMAP GENERATION (HEAT)
    ###################################################################
    ###################################################################

    ## HEAT - header (2) - Heatmap analysis
    output$headheat <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Interactive Heatmap (click on cells)")
        }
    })

    ## HEAT - input - Choose number of variable IDs
    output$heatnumber <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            textInput(
                inputId = "heatnumber",
                label = "ID cutoff",
                value = 20
            )
        }
    })

    ## HEAT - data - get variable IDs
    heattran2 <- reactive({
        if (input$goqc == 0) {
            return()
        } else {
            dds.counts <- ddsout()[[3]]
            heat.counts <- ddstran()[[1]]
            heat.counts <- assay(heat.counts)
            num <- input$heatnumber
            topID <- order(rowMeans(dds.counts), decreasing = TRUE)
            heat.mat <- heat.counts[topID, ]
            # heat.mat <- heat.mat - rowMeans(heat.mat)
            heat.mat <- heat.mat[1:num, ,drop = FALSE]
        }
        return(list(heat.mat))
    })

    ## HEAT - visualization - Plotly heatmap
    output$heatplot1 <- renderPlotly({
        validate(
            need(input$heatnumber != "", "")
        )
        if (input$goqc == 0) {
            return()
        } else {
            heat <- heattran2()[[1]]
            tooltips <- paste0(
                # "<b>Sample:</b> ", colnames(heat), "<br />",
                "<b>ID:</b> ", rownames(heat), "<br />",
                "<b>Count:</b> ", round(heat, 3)
            )
            tooltips <- matrix(tooltips, ncol = ncol(heat), byrow = FALSE)
            plot_ly(
                x = colnames(heat),
                y = rownames(heat),
                z = heat,
                type = "heatmap",
                text = tooltips,
                hoverinfo = "text",
                source = "heatplot"
            ) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = ""),
                margin = list(l = 100)
            )
        }
    })

    ## HEAT - Show download button - Heatmap (PDF)
    output$dlqcheatplot1pdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcheatplot1pdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## HEAT - Download plot - Heatmap (PDF)
    output$dlqcheatplot1pdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-heatmap.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 7.5, onefile = FALSE) # open the pdf device
            qcHeatMap(
                heat = heattran2()[[1]],
                n = input$heatnumber
            )
            dev.off()
        }
    )

    ## HEAT - Show download button - Heatmap (PNG)
    output$dlqcheatplot1png <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcheatplot1pngimg", "Download Static R Plot (PNG)")
        }
    })

    ## HEAT - Download plot - Heatmap (PNG)
    output$dlqcheatplot1pngimg <- downloadHandler(
        filename =  function() {
            paste("qc-heatmap.png")
        },
        content = function(file) {
            png(file, width = 800, height = 850)
            qcHeatMap(
                heat = heattran2()[[1]],
                n = input$heatnumber
            )
            dev.off()
        }
    )

    ##HEAT - select input - choose factor - HEATMAP
    output$heatfactor <- renderUI({
        s <- event_data("plotly_click", source = "heatplot")
        validate(
            need(
                s != "",
                message = "")
        )
        tmp <- ddsout()[[2]]
        selectInput(
            inputId = "heatfactor",
            label = "Choose factor",
            choices = colnames(tmp)
        )
    })

    ## HEAT - visualization - count plot
    output$heatplot2 <- renderPlotly({
        if (input$goqc == 0) {
            return()
        } else {
            s <- event_data("plotly_click", source = "heatplot")
            validate(
                need(
                    s != "",
                    message = "Click on one of the heatmap cells to view this plot!"
                )
            )
            rc.data <- counts(ddsout()[[1]])
            test <- getGenes(
                rc.data = rc.data,
                id = s[["y"]],
                coldata = ddsout()[[2]]
            )
            tooltips <- paste0(
                "<b>Sample:</b> ", test$sample, "<br />",
                "<b>Counts:</b> ", round(test$counts, 3)
            )

            plot_ly(
                data = test,
                type = "scatter",
                mode = "markers",
                x = test[, input$heatfactor],
                y = test[, "counts"],
                color = test[, input$heatfactor],
                text = tooltips,
                marker = list(size = 7),
                hoverinfo = "text"
            ) %>%
            layout(
                title = paste(s[["y"]], "Counts"),
                xaxis = list(title = paste(input$fact)),
                yaxis = list(title = "Normalized counts")
            )
        }
    })

    ## HEAT - Show download button - Heat counts (PDF)
    output$dlqcheatplot2pdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcheatplot2pdfimg", "Download Static R Plot (PDF)")
        }
    })

    ## HEAT - Download plot - Heat counts (PDF)
    output$dlqcheatplot2pdfimg <- downloadHandler(
        filename =  function() {
            paste("qc-heat-counts.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
            qcHeatCount(
                s = event_data("plotly_click", source = "heatplot"),
                rc.data = counts(ddsout()[[1]]),
                coldata = ddsout()[[2]],
                heatfactor = input$heatfactor
            )
            dev.off()
        }
    )

    ## HEAT - Show download button - Heat counts (PNG)
    output$dlqcheatplot2png <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcheatplot2pngimg", "Download Static R Plot (PNG)")
        }
    })

    ## HEAT - Download plot - Heat counts (PNG)
    output$dlqcheatplot2pngimg <- downloadHandler(
        filename =  function() {
            paste("qc-heat-counts.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            qcHeatCount(
                s = event_data("plotly_click", source = "heatplot"),
                rc.data = counts(ddsout()[[1]]),
                coldata = ddsout()[[2]],
                heatfactor = input$heatfactor
            )
            dev.off()
        }
    )





    ###################################################################
    ###################################################################
    ### SECTION TMP - CLUSTERING ALGORITHMS (CLUST)
    ###################################################################
    ###################################################################

    ## CLUST - header (2) - clustering
    output$headclust <- renderUI({
        if(input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Clustering Analysis")
        }
    })

    output$headclustwarn <- renderUI({
        if (input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            p(
                br(),
                em(
                    paste0(
                        "WARNING: running these algorithms may be very ",
                        "slow, depending on the number of genes in analysis. ",
                        "In case of server timeout, consider ",
                        "downloading the IRIS source and run locally."
                    )
                )
            )
        }
    })

    ## BIC - input - Choose number of variable genes
    output$clustvarnumber <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            textInput(
                inputId = "clustvarnumber",
                label = "Variable cutoff value",
                value = 500
            )
        }
    })

    ## CLUST - input - Choose clustering algorithm
    output$clustalg <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            selectInput(
                inputId = "clustalg",
                label = "Choose clustering algorithm",
                choices = c(
                    "WGCNA" = "wgcna",
                    "K-Medoids" = "kmed",
                    "MCL" = "mcl"
                ),
                selected = "WGCNA"
            )
        }
    })

    ## CLUST - actionbutton - submit clustering
    output$goclust <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            actionButton(
                "goclust",
                "Launch Analysis",
                icon = icon("space-shuttle")
            )
        }
    })

    ## CLUST - reactive - get variable counts (WIP - very mess atm...)
    clustout <- eventReactive(input$goclust, {
        num <- input$clustvarnumber
        cts <- ddsout()[[1]]
        cts <- assay(cts)
        tran <- ddstran()[[1]]
        tran <- assay(tran)
        topID <- order(rowVars(cts), decreasing = TRUE)
        cts.var <- tran[topID, ]
        dds_mat <- cts.var[1:num, ]
        if (input$clustalg == "wgcna") {
            withProgress(message = "Running WGCNA...", value = 0, {
                incProgress(1/4)
                enableWGCNAThreads()

                dds_mat <- as.matrix(dds_mat)
                gene.names <- sort(rownames(dds_mat))

                datExpr <- t(log2(dds_mat + 1))

                incProgress(2/4)
                # Run this to check if there are gene outliers
                gsg <- goodSamplesGenes(datExpr, verbose = 3)
                gsg$allOK

                # Create an object called "datTraits" that contains your
                # trait data
                datTraits <- colData(ddsout()[[1]])
                head(datTraits)

                # Form a data frame analogous to expression data that will
                # hold the clinical traits.
                # should return TRUE if datasets align correctly, otherwise your
                # names are out of order
                table(rownames(datTraits) == rownames(datExpr))
                A <- adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
                k <- as.numeric(apply(A,2,sum))-1 # standardized connectivity
                Z.k <- scale(k)
                thresholdZ.k <- -2.5 # often -2.5
                outlierColor <- ifelse(Z.k<thresholdZ.k,"red","black")
                sampleTree <- flashClust(as.dist(1-A), method = "average")

                # Convert traits to a color representation where red indicates high values
                traitColors <- data.frame(labels2colors(datTraits))
                dimnames(traitColors)[[2]] <- paste(names(datTraits))
                datColors <- data.frame(outlier = outlierColor, traitColors)

                incProgress(3/4)
                # TOM analysis - (computationally expensive)
                enableWGCNAThreads()
                softPower <- 18
                adjacency <- adjacency(datExpr, power = softPower, type = "signed") #specify network type
                TOM <- TOMsimilarity(adjacency, TOMType = "signed") # specify network type
                dissTOM <- 1 - TOM
                geneTree <- flashClust(as.dist(dissTOM), method="average")

                # This sets the minimum number of genes to cluster into a module
                minModuleSize = 30
                dynamicMods <- cutreeDynamic(
                    dendro = geneTree,
                    distM = dissTOM,
                    deepSplit = 2,
                    pamRespectsDendro = FALSE,
                    minClusterSize = minModuleSize
                )

                dynamicColors <- labels2colors(dynamicMods)
                MEList <- moduleEigengenes(
                    datExpr,
                    colors = dynamicColors,
                    softPower = 18
                )
                MEs <- MEList$eigengenes
                MEDiss <- 1 - cor(MEs)
                METree <- flashClust(as.dist(MEDiss), method = "average")

                # set a threhold for merging modules. In this example we are
                # not merging so MEDissThres=0.0
                MEDissThres <- 0.0
                merge <- mergeCloseModules(
                    datExpr,
                    dynamicColors,
                    cutHeight = MEDissThres,
                    verbose = 3
                )

                mergedColors <- merge$colors
                mergedMEs <- merge$newMEs

                # Set the diagonal of the dissimilarity to NA
                diag(dissTOM) = NA;

                # Export modules to data frame
                module_colors= setdiff(unique(dynamicColors), "grey")
                modlist <- list()
                for (color in module_colors){
                    module = gene.names[which(dynamicColors == color)]
                    modlist[[color]] <- list(
                        gene = module
                    )
                }

                moddf <- data.frame(unlist(modlist))
                moddf$module <- gsub("\\..*", "", row.names(moddf))
                colnames(moddf)[1] <- "gene"
                rownames(moddf) <- seq_len(nrow(moddf))
                moddf$gene <- as.character(moddf$gene)
                moddf$module <- as.factor(moddf$module)

                sampleDF <- data.frame(
                    sample = sampleTree$labels,
                    outlier = datColors$outlier,
                    cluster = datColors$Cluster
                )

                return(
                    list(
                        sampleTree,
                        datColors,
                        geneTree,
                        dynamicColors,
                        mergedColors,
                        dissTOM,
                        moddf,
                        sampleDF
                    )
                )
                incProgress(4/4)
            })
        } else if (input$clustalg == "kmed") {
            withProgress(message = "Running K-Medoids...", value = 0, {
                incProgress(1/2)
                num <- as.matrix(dds_mat)
                mrwdist <- distNumeric(num, num, method = "mrw")

                # Detect cores of machine
                nclust <- parallel::detectCores() - 1
                message("Using ", nclust, " threads...")
                result <- fastkmed(mrwdist, ncluster = nclust, iterate = 50)

                # a simple and fast k-medoids function for bootstrap evaluation
                parkboot <- function(x, nclust) {
                    res <- fastkmed(x, nclust, iterate = 50)
                    return(res$cluster)
                }

                fastkmedboot <- clustboot(mrwdist, nclust = nclust, parkboot, nboot = 50)

                # consensus matrix
                wardorder <- function(x, nclust) {
                    res <- fastcluster::hclust(x, method = "ward.D2")
                    member <- cutree(res, nclust)
                    return(member)
                }
                consensusfastkmed <- consensusmatrix(fastkmedboot, nclust = nclust, wardorder)

                # data frame generation
                output <- data.frame(
                    gene_id = rownames(num),
                    cluster = result$cluster
                )
                rownames(output) <- seq_len(nrow(output))
                output$gene_id <- as.character(output$gene_id)
                output$cluster <- as.factor(output$cluster)

                return(
                    list(
                        consensusfastkmed,
                        output
                    )
                )
                incProgress(2/2)
            })
        } else if (input$clustalg == "mcl") {
            withProgress(message = "Running MCL...", value = 0, {
                incProgress(1/2)
                message("Running MCL...")
                exp_file <- dds_mat
                row_sub = apply(exp_file, 1, function(row) all(row != 0))

                # data pre-processing replace NAs with row means
                for (i in 1:ncol(exp_file)) {
                    exp_file[is.na(exp_file[, i]), i] <- mean(exp_file[, i], na.rm = TRUE)
                }
                exp_file <- exp_file[apply(exp_file, 1, function(x) !all(x == 0)), ]  # remove rows that are all 0
                exp_cor <- cor(t(exp_file), method = "spearman", use = "complete.obs")  # gene to gene correlation matrix
                MCL_name <- rownames(exp_cor)  # the vertix
                Covered <- length(MCL_name)  # the #of covered cells

                calculate_mcl <- function(i) {
                    # main mcl function, used in apply function to speed up calculation
                    cluster <- list()
                    cluster <- mcl(exp_cor, addLoops = F, inflation = 50, max.iter = 500)
                    return(cluster)
                }

                i <- as.data.frame(seq(5))  # test mcl inflation parameter from 1 to 50
                i <- as.data.frame(50)  # may be just set inflation=50? one iteration is about 5 minutes

                cluster_all <- apply(i, 1, calculate_mcl)


                KK <- as.data.frame(do.call(rbind, lapply(cluster_all, "[[", 1)))  # extract the number of clusters
                CAN_I <- c(which(as.numeric(as.character(KK$V1)) >= 2))  # results that has more than 5 clusters
                tt <- as.numeric(as.character(KK$V1))
                tt <- sort(table(tt), decreasing = T)[1]
                Final_K <- as.numeric(names(tt))

                if (Final_K == 1) {
                    message("Final K is 1...")
                    matdist <- 0
                    result <- data.frame()
                    hc <- 0
                } else {
                    if (length(CAN_I) != 0) {
                        message("Final K is NOT 1 - continuing...")
                        MATRIX <- rep(0, Covered) %o% rep(0, Covered)
                        for (k in 1:length(CAN_I)) {
                            MCL_label <- cluster_all[[CAN_I[k]]]$Cluster  # record the label
                            ClusterNum <- unique(MCL_label)  # record the number of clusters
                            TEMP <- rep(0, Covered) %o% rep(0, Covered)
                            temp <- rep(0, Covered) %o% rep(0, length(ClusterNum))
                            for (n in 1:length(ClusterNum)) {
                                index <- which(MCL_label == ClusterNum[n])
                                temp[index, n] <- 1
                                TEMP <- TEMP + temp[, n] %o% temp[, n]
                            }
                            MATRIX <- MATRIX + TEMP
                        }
                        MATRIX <- MATRIX/length(CAN_I)
                        rownames(MATRIX) <- colnames(MATRIX) <- rownames(exp_cor)
                        hc <- hclust(dist(MATRIX))
                        memb <- cutree(hc, k = Final_K)
                    }
                    if (length(rownames(exp_cor)) == length(rownames(exp_cor))) {
                        label <- memb
                    } else {
                        LEFT <- setdiff(names(RAW)[-1], MCL_name)
                        LEFT_Cluster <- rep(Final_K + 1, length(LEFT))
                        df_cell_label <- data.frame(
                            cell = c(names(memb), LEFT),
                            cluster = c(
                                memb,
                                LEFT_Cluster
                            ),
                            K = rep(Final_K + 1, length(rownames(exp_cor)))
                        )
                        label <- df_cell_label$cluster
                    }
                    matdist <- dist(MATRIX)
                    matdist <- data.matrix(matdist)
                    result <- as.data.frame(label)
                    result$cluster <- result$label
                    result$gene_id <- rownames(result)
                    result <- result[, c(3, 2)]
                    rownames(result) <- seq_len(nrow(result))
                }


                message("Sending data to list...")
                return(
                    list(
                        matdist,
                        result,
                        hc,
                        Final_K
                    )
                )
                incProgress(2/2)
            })
        }
    })

    ## CLUST - head (1) - sample dendrogram
    output$headclustplotW01 <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            h5(strong("WGCNA - Sample Dendrogram"))
        } else if (input$clustalg == "kmed") {
            h5(strong("K-Medoids - Consensus Matrix Heatmap"))
        } else if (input$clustalg == "mcl") {
            Final_K <- clustout()[[4]]
            if (Final_K == 1) {
                h5(strong("MCL - Cluster Diagram - No Clusters Found!"))
            } else {
                h5(strong("MCL - Cluster Diagram"))
            }
        }
    })

    ## CLUST - visualize (WGCNA) - dendrogram plots 1
    output$clustplotW01 <- renderPlot({
        req(clustout())
        if (input$clustalg == "wgcna") {
            validate(
                need(input$goqc != "", "")
            )
            sampleTree <- clustout()[[1]]
            datColors <- clustout()[[2]]

            plotDendroAndColors(
                sampleTree,
                groupLabels = names(datColors),
                colors = datColors,
                main = "Sample Dendrogram and Module Colors"
            )
        } else if (input$clustalg == "kmed") {
            validate(
                need(input$goqc != "", "")
            )
            consensusfastkmed <- clustout()[[1]]
            clustheatmap(
                consensusfastkmed,
                "K-Medoids Consensus Matrix Heatmap"
            )
        } else if (input$clustalg == "mcl") {
            validate(
                need(input$goqc != "", "")
            )
            matdist <- clustout()[[1]]
            hc <- clustout()[[3]]
            Final_K <- clustout()[[4]]
            if (Final_K == 1) {
                return()
            } else {
                plot(as.phylo(hc), type = "fan", cex = 0.3)
            }
        }
    })

    ## CLUST - sample dendrogram download (PNG) 1
    output$downloadclustplotW01png <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW01pngimg",
                "Download Plot (PNG)"
            )
        } else if (input$clustalg == "kmed") {
            downloadButton(
                "downloadclustplotK01pngimg",
                "Download Plot (PNG)"
            )
        } else if (input$clustalg == "mcl") {
            Final_K <- clustout()[[4]]
            if (Final_K == 1) {
                return()
            } else {
                downloadButton(
                    "downloadclustplotM01pngimg",
                    "Download Plot (PNG)"
                )
            }
        }
    })

    ## CLUST - sample dendrogram download (PNG) 2
    output$downloadclustplotW01pngimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-sample-dendrogram.png")
        },
        content = function(file) {
            png(file, width = 800, height = 600)
            sampleTree <- clustout()[[1]]
            datColors <- clustout()[[2]]

            plotDendroAndColors(
                sampleTree,
                groupLabels = names(datColors),
                colors = datColors,
                main = "Sample Dendrogram and Module Colors"
            )
            dev.off()
        }
    )

    ## CLUST - download (2) - consensus matrix
    output$downloadclustplotK01pngimg <- downloadHandler(
        filename = function() {
            paste("clust-kmed-consensus-matrix.png")
        },
        content = function(file) {
            req(clustout())
            png(file, width = 800, height = 600)
            sample <- clustout()[[1]]
            output <- clustheatmap(sample, "K-Medoids Consensus Matrix Heatmap")
            print(output)
            dev.off()
        }
    )

    ## CLUST - download (3) - distance matrix (MCL)
    output$downloadclustplotM01pngimg <- downloadHandler(
        filename = function() {
            paste("clust-mcl-cluster-diagram.png")
        },
        content = function(file) {
            req(clustout())
            matdist <- clustout()[[1]]
            hc <- clustout()[[3]]
            png(file, width = 800, height = 600)
            output <- plot(as.phylo(hc), type = "fan", cex = 0.2)
            print(output)
            dev.off()
        }
    )

    ## CLUST - sample dendrogram download (pdf) 1
    output$downloadclustplotW01pdf <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW01pdfimg",
                "Download Plot (PDF)"
            )
        } else if (input$clustalg == "kmed") {
            downloadButton(
                "downloadclustplotK01pdfimg",
                "Download Plot (PDF)"
            )
        } else if (input$clustalg == "mcl") {
            Final_K <- clustout()[[4]]
            if (Final_K == 1) {
                return()
            } else {
                downloadButton(
                    "downloadclustplotM01pdfimg",
                    "Download Plot (PDF)"
                )
            }
        }
    })

    ## CLUST - sample dendrogram download (pdf) 2
    output$downloadclustplotW01pdfimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-sample-dendrogram.pdf")
        },
        content = function(file) {
            pdf(file, width = 12, height = 8)
            sampleTree <- clustout()[[1]]
            datColors <- clustout()[[2]]

            plotDendroAndColors(
                sampleTree,
                groupLabels = names(datColors),
                colors = datColors,
                main = "Sample Dendrogram and Module Colors"
            )
            dev.off()
        }
    )

    ## CLUST - download (3) - distance matrix (MCL) - pdf
    output$downloadclustplotM01pdfimg <- downloadHandler(
        filename = function() {
            paste("clust-mcl-cluster-diagram.pdf")
        },
        content = function(file) {
            req(clustout())
            matdist <- clustout()[[1]]
            hc <- clustout()[[3]]
            pdf(file, width = 8, height = 6)
            plot(as.phylo(hc), type = "fan", cex = 0.2)
            dev.off()
        }
    )

    ## CLUST - head (2) - gene dendrogram
    output$headclustplotW02 <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            h5(strong("WGCNA - Gene Dendrogram"))
        } else {
            return()
        }
    })

    ## CLUST - visualize (WGCNA) - gene dendrogram
    output$clustplotW02 <- renderPlot({
        req(clustout())
        if (input$clustalg == "wgcna") {
            validate(
                need(input$goqc != "", "")
            )
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            mergedColors <- clustout()[[5]]

            plotDendroAndColors(
                geneTree,
                cbind(dynamicColors, mergedColors),
                c("Dynamic Tree Cut", "Merged dynamic"),
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05,
                main = "Gene Dendrogram and Module Colors"
            )
        } else {
            return()
        }
    })

    ## CLUST - gene dendrogram download (PNG) 1
    output$downloadclustplotW02png <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW02pngimg",
                "Download Plot (PNG)"
            )
        } else {
            return()
        }
    })

    ## CLUST - gene dendrogram download (PNG) 2
    output$downloadclustplotW02pngimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-gene-dendrogram.png")
        },
        content = function(file) {
            png(file, width = 1000, height = 800)
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            mergedColors <- clustout()[[5]]

            plotDendroAndColors(
                geneTree,
                cbind(dynamicColors, mergedColors),
                c("Dynamic Tree Cut", "Merged dynamic"),
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05,
                main = "Gene Dendrogram and Module Colors"
            )
            dev.off()
        }
    )

    ## CLUST - gene dendrogram download (pdf) 1
    output$downloadclustplotW02pdf <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW02pdfimg",
                "Download Plot (pdf)"
            )
        } else {
            return()
        }
    })

    ## CLUST - gene dendrogram download (pdf) 2
    output$downloadclustplotW02pdfimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-gene-dendrogram.pdf")
        },
        content = function(file) {
            pdf(file, width = 12, height = 8)
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            mergedColors <- clustout()[[5]]

            plotDendroAndColors(
                geneTree,
                cbind(dynamicColors, mergedColors),
                c("Dynamic Tree Cut", "Merged dynamic"),
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05,
                main = "Gene Dendrogram and Module Colors"
            )
            dev.off()
        }
    )

    ## CLUST - head (3) - sample dendrogram
    output$headclustplotW03 <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            h5(strong("WGCNA - Topological Overlap Matrix"))
        } else {
            return()
        }
    })

    ## CLUST - visualize (WGCNA) - TOM plot
    output$clustplotW03 <- renderPlot({
        req(clustout())
        if (input$clustalg == "wgcna") {
            validate(
                need(input$goqc != "", "")
            )
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            dissTOM <- clustout()[[6]]

            TOMplot(
                dissTOM^4,
                geneTree,
                as.character(dynamicColors),
                main = "TOM Plot"
            )
        } else {
            return()
        }
    })

    ## CLUST - TOM plot download (PNG) 1
    output$downloadclustplotW03png <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW03pngimg",
                "Download Plot (PNG)"
            )
        } else {
            return()
        }
    })

    ## CLUST - TOM plot download (PNG) 2
    output$downloadclustplotW03pngimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-tom-plot.png")
        },
        content = function(file) {
            png(file, width = 800, height = 600)
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            dissTOM <- clustout()[[6]]

            TOMplot(
                dissTOM^4,
                geneTree,
                as.character(dynamicColors),
                main = "TOM Plot"
            )
            dev.off()
        }
    )

    ## CLUST - TOM plot download (pdf) 1
    output$downloadclustplotW03pdf <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustplotW03pdfimg",
                "Download Plot (pdf)"
            )
        } else {
            return()
        }
    })

    ## CLUST - TOM plot download (pdf) 2
    output$downloadclustplotW03pdfimg <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-tom-plot.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6)
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            dissTOM <- clustout()[[6]]

            TOMplot(
                dissTOM^4,
                geneTree,
                as.character(dynamicColors),
                main = "TOM Plot"
            )
            dev.off()
        }
    )

    ## CLUST - head (4) - download gene modules
    output$headclustmoddown <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            h5(strong("WGCNA - Download Gene Modules"))
        } else {
            return()
        }
    })

    ## CLUST - gene module download 1
    output$downloadclustmod <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustmod2",
                "Download Gene Modules (CSV)"
            )
        } else {
            return()
        }
    })

    ## CLUST - gene module download 2
    output$downloadclustmod2 <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-gene-modules.csv")
        },
        content = function(file) {
            moddf <- clustout()[[7]]
            write.csv(moddf, file, row.names = FALSE, col.names = TRUE)
        }
    )

    ## CLUST - sample module download 1
    output$downloadclustsample <- renderUI({
        req(clustout())
        if (input$clustalg == "wgcna") {
            downloadButton(
                "downloadclustsample2",
                "Download Sample Modules (CSV)"
            )
        } else {
            return()
        }
    })

    ## CLUST - sample module download 2
    output$downloadclustsample2 <- downloadHandler(
        filename = function() {
            paste("clust-wgcna-sample-modules.csv")
        },
        content = function(file) {
            sampledf <- clustout()[[8]]
            write.csv(sampledf, file, row.names = FALSE, col.names = TRUE)
        }
    )

    ## CLUST - download (4) - consensus matrix
    output$downloadclustplotK01pdfimg <- downloadHandler(
        filename = function() {
            paste("clust-kmed-consensus-matrix.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6)
            consensusfastkmed <- clustout()[[1]]
            p <- clustheatmap(
                consensusfastkmed,
                "K-Medoids Consensus Matrix Heatmap"
            )
            print(p)
            dev.off()
        }
    )

    ## CLUST - gene module download 1
    output$downloadclustmodK <- renderUI({
        req(clustout())
        if (input$clustalg == "kmed") {
            downloadButton(
                "downloadclustmodK2",
                "Download Clusters (CSV)"
            )
        } else if (input$clustalg == "mcl") {
            Final_K <- clustout()[[4]]
            if (Final_K == 1) {
                return()
            } else {
                downloadButton(
                    "downloadclustmodM2",
                    "Download Clusters (CSV)"
                )
            }
        } else {
            return()
        }
    })

    ## CLUST - gene module download 2
    output$downloadclustmodK2 <- downloadHandler(
        filename = function() {
            paste("clust-kmed-gene-clusters.csv")
        },
        content = function(file) {
            clustdf <- clustout()[[2]]
            write.csv(clustdf, file, row.names = FALSE, col.names = TRUE)
        }
    )

    ## CLUST - gene module download 2
    output$downloadclustmodM2 <- downloadHandler(
        filename = function() {
            paste("clust-mcl-gene-clusters.csv")
        },
        content = function(file) {
            clustdf <- clustout()[[2]]
            write.csv(clustdf, file, row.names = FALSE, col.names = TRUE)
        }
    )

    ###################################################################
    ###################################################################
    ### SECTION 05 - BICLUSTERING ALGORITHMS (BIC)
    ###################################################################
    ###################################################################

    ## BIC - header (2) - Bicluster analysis
    output$headbic <- renderUI({
        if (input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Bicluster Analysis")
        }
    })

    ## BIC - header (3) - Bicluster analysis parameters
    output$headbicparameters <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            h5(strong("Parameters"))
        }
    })

    ## BIC - input - Choose number of variable genes
    output$bicvarnumber <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            textInput(
                inputId = "bicvarnumber",
                label = "Variable cutoff value",
                value = 500
            )
        }
    })

    ## BIC - input - Choose bicluster algorithm
    output$bicalg <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            selectInput(
                inputId = "bicalg",
                label = "Choose bicluster algorithm",
                choices = c(
                    "QUBIC" = "qubic",
                    "Bimax" = "bimax",
                    "CC" = "cc",
                    "Plaid" = "plaid",
                    "Spectral" = "spectral",
                    "Xmotifs" = "xmotifs"
                )
            )
        }
    })

    ## BIC - actionbutton - submit biclustering
    output$gobic <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            actionButton(
                "gobic",
                "Launch Analysis",
                icon = icon("space-shuttle")
            )
        }
    })

    ## BIC - reactive - get variable counts
    bicout <- eventReactive(input$gobic, {
        num <- input$bicvarnumber
        cts <- ddsout()[[1]]
        cts <- assay(cts)
        tran <- ddstran()[[1]]
        tran <- assay(tran)
        topID <- order(rowVars(cts), decreasing = TRUE)
        cts.var <- tran[topID, ]
        cts.var <- cts.var[1:num, ]
        if (input$bicalg == "qubic") {
            withProgress(message = "Running QUBIC...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCQU() )
                incProgress(2/2)
            })
        } else if (input$bicalg == "bimax") {
            withProgress(message = "Running Bimax...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCBimax() )
                incProgress(2/2)
            })
        } else if (input$bicalg == "cc") {
            withProgress(message = "Running CC...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCCC() )
                incProgress(2/2)
            })
        } else if (input$bicalg == "plaid") {
            withProgress(message = "Running Plaid...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCPlaid() )
                incProgress(2/2)
            })
        } else if (input$bicalg == "spectral") {
            withProgress(message = "Running Spectal...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCSpectral() )
                incProgress(2/2)
            })
        } else if (input$bicalg == "xmotifs") {
            withProgress(message = "Running Xmotifs...", value = 0, {
                incProgress(1/2)
                res <- biclust::biclust(x = cts.var, method = BCXmotifs() )
                incProgress(2/2)
            })
        }
        return(list(res, cts.var))
    })

    ## BIC - input - choose cluster
    output$bicclustnumber <- renderUI({
        res <- bicout()[[1]]
        if (res@Number < 1) {
            return()
        } else {
            n <- res@Number
            selectInput(
                inputId = "bicclustnumber",
                label = "Choose cluster",
                choices = c(1:n)
            )
        }
    })

    ## BIC - head (3) - Clust overview
    output$headbicsummary <- renderUI({
        validate(
            need(
                input$bicclustnumber != "",
                ""
            )
        )
        h5(strong("Clustering Overview"))
    })


    ## BIC - text - cluster overview
    output$bicsummary <- renderUI({
        res <- bicout()[[1]]
        if (res@Number < 1) {
            p("No clusters found using this algorithm!")
        } else {
            rn <- res@RowxNumber
            cn <- res@NumberxCol
            clust <- input$bicclustnumber
            clust <- as.numeric(clust)
            ids <- length(rn[, clust][rn[, clust] == TRUE])
            samp <- length(cn[clust, ][cn[clust, ] == TRUE])

            out <- paste0(
                "This algorithm found ", ids, " IDs amongst ", samp,
                " samples in cluster ", clust,
                ". To see what IDs were found in this cluster, click on the ",
                "'Download Cluster IDs' button at the bottom of the page."
            )
            p(paste(out))
        }
    })

    ## BIC - head (3) - Clust overview
    output$headbicsummary <- renderUI({
        validate(
            need(
                input$bicclustnumber != "",
                ""
            )
        )
        h5(strong("Clustering Overview"))
    })

    ## BIC - header (3) - Bicluster analysis parameters
    output$headbicheatsummary <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            res <- bicout()[[1]]
            if (res@Number < 1) {
                return()
            } else {
                clust <- input$bicclustnumber
                h5(strong(paste0("Heatmap analysis of cluster ", clust)))
            }
        }
    })

    ## BIC - visualize - Bicluster heatmap
    output$bicheatplot <- renderPlot({
        validate(
            need(input$bicclustnumber != "", "")
        )
        res <- bicout()[[1]]
        if (res@Number < 1) {
            return()
        } else {
            bicPlot(
                n = input$bicclustnumber,
                res = bicout()[[1]],
                cts.var = bicout()[[2]]
            )
        }
    })

    ## BIC - reactive - get cluster names
    bicout2 <- eventReactive(input$gobic, {
        res <- bicout()[[1]]
        cts.var <- bicout()[[2]]
        rn <- res@RowxNumber
        col.num <- as.numeric(input$bicclustnumber)
        select <- rn[, col.num]
        cts.nam <- cts.var[which(select), ]
        cts.nam <- row.names(cts.nam)
        return(list(cts.nam))
    })

    ## BIC - Show download button (FILTERED)
    output$downloadbicfilt <- renderUI({
        validate(
            need(
                expr = !is.null(input$gobic),
                message = ""
            )
        )
        if(input$gobic == 0) {
            return()
        } else {
            res <- bicout()[[1]]
            if (res@Number < 1) {
                return()
            } else {
                downloadButton("downloadbicclust", "Download Cluster IDs")
            }
        }
    })

    ## BIC - Download data
    output$downloadbicclust <- downloadHandler(
        filename = function() {
            n <- as.numeric(input$bicclustnumber)
            paste0("cluster-", n, "-ids.csv")
        },
        content = function(file) {
            cts.nam <- bicout2()[[1]]
            cts.nam <- data.frame(IDs = cts.nam)
            write.csv(cts.nam, file, row.names = FALSE, col.names = TRUE)
        }
    )

    ## BIC - Show download button (FILTERED)
    output$downloadbicplotpdf <- renderUI({
        validate(
            need(
                expr = !is.null(input$gobic),
                message = ""
            )
        )
        if(input$gobic == 0) {
            return()
        } else {
            res <- bicout()[[1]]
            if (res@Number < 1) {
                return()
            } else {
                downloadButton("downloadbicplotpdfimg", "Download Plot (PDF)")
            }
        }
    })

    ## BIC - Download plot
    output$downloadbicplotpdfimg <- downloadHandler(
        filename =  function() {
            paste("biclust-plot.pdf")
        },
        content = function(file) {
            pdf(file, width = 6, height = 8) # open the pdf device
            bicPlot(
                n = input$bicclustnumber,
                res = bicout()[[1]],
                cts.var = bicout()[[2]]
            )
            dev.off()
        }
    )

    ## BIC - Show download button (FILTERED) (PNG)
    output$downloadbicplotpng <- renderUI({
        validate(
            need(
                expr = !is.null(input$gobic),
                message = ""
            )
        )
        if(input$gobic == 0) {
            return()
        } else {
            res <- bicout()[[1]]
            if (res@Number < 1) {
                return()
            } else {
                downloadButton("downloadbicplotpngimg", "Download Plot (PNG)")
            }
        }
    })

    ## BIC - Download plot (PNG)
    output$downloadbicplotpngimg <- downloadHandler(
        filename =  function() {
            paste("biclust-plot.png")
        },
        content = function(file) {
            png(file, width = 600, height = 800)
            bicPlot(
                n = input$bicclustnumber,
                res = bicout()[[1]],
                cts.var = bicout()[[2]]
            )
            dev.off()
        }
    )





    ###################################################################
    ###################################################################
    ### SECTION 06 - CORRELATION ANALYSIS (COR)
    ###################################################################
    ###################################################################

    ## COR - header (2) - Correlation analysis
    output$headcor <- renderUI({
        if (input$goqc == 0) {
            p(
                br(),
                em(
                    "Load data and click the 'submit' button on the 'Submit and QC' tab to see the results."
                ),
                style = "color:grey"
            )
        } else {
            h4("Interactive Correlation Analysis (click on cells)")
        }
    })

    ## COR - reactive - correlation matrix
    corout <- eventReactive(input$goqc, {
        cts <- ddstran()[[1]]
        cts <- assay(cts)
        cor.mat <- cor(cts)
        return(list(cor.mat, cts))
    })


    ## COR - visaulization - correlation matrix (Plotly)
    output$corplot1 <- renderPlotly({
        validate(
            need(
                expr = !is.null(corout()),
                message = ""
            )
        )
        withProgress(message = "Rendering correlation matrix...", value = 0, {
            if (input$goqc == 0) {
                return()
            } else {
                incProgress()
                cor.mat <- corout()[[1]]
                cor.mat <- as.matrix(cor.mat)
                tooltips <- paste0(
                    "<b>R:</b> ", round(cor.mat, 3)
                )
                tooltips <- matrix(
                    tooltips,
                    ncol = ncol(cor.mat),
                    byrow = FALSE
                )
                plot_ly(
                    x = colnames(cor.mat),
                    y = rownames(cor.mat),
                    z = cor.mat,
                    type = "heatmap",
                    text = tooltips,
                    hoverinfo = "text",
                    source = "corplot"
                ) %>%
                layout(
                    xaxis = list(title = ""),
                    yaxis = list(title = ""),
                    margin = list(l = 100)
                )
            }
        })
    })

    ## COR - DOWNLOAD BUTTON - Corplot1 (pdf)
    output$dlqccorplot1pdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                "dlqccorplot1pdfimg",
                "Download Static R Plot (PDF)"
            )
        }
    })

    ## COR - DOWNLOAD PLOT - Corplot1 (pdf)
    output$dlqccorplot1pdfimg <- downloadHandler(
        filename =  function() {
            paste("cor-matrix.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 6.5, onefile = FALSE)
            corMatPlot(
                cor.mat = corout()[[1]]
            )
            dev.off()
        }
    )

    ## COR - DOWNLOAD BUTTON - Corplot1 (png)
    output$dlqccorplot1png <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                "dlqccorplot1pngimg",
                "Download Static R Plot (PNG)"
            )
        }
    })

    ## COR - DOWNLOAD PLOT - Corplot1 (png)
    output$dlqccorplot1pngimg <- downloadHandler(
        filename =  function() {
            paste("cor-matrix.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            corMatPlot(
                cor.mat = corout()[[1]]
            )
            dev.off()
        }
    )

    ## COR - visualization - scatterplot
    output$corplot2 <- renderPlotly({
        if (input$goqc == 0) {
            return()
        } else {
            withProgress(message = "Rendering count plot...", value = 0, {
                incProgress()
                s.cor <- event_data("plotly_click", source = "corplot")
                validate(
                    need(
                        s.cor != "",
                        message = "Click on one of the heatmap cells to view this plot!"
                    )
                )
                cts.tran <- corout()[[2]]
                cts.tran <- as.data.frame(cts.tran)
                x <- s.cor[["x"]]
                y <- s.cor[["y"]]
                tooltips <- paste0(
                    "<b>ID:</b> ", rownames(cts.tran)
                )
                plot_ly(
                    data = cts.tran,
                    type = "scatter",
                    mode = "markers",
                    x = cts.tran[, x],
                    y = cts.tran[, y],
                    text = tooltips,
                    marker = list(size = 2),
                    hoverinfo = "text"
                ) %>%
                layout(
                    title = paste(y, "vs.", x),
                    xaxis = list(title = x),
                    yaxis = list(title = y)
                )
            })
        }
    })

    ## COR - DOWNLOAD BUTTON - corplot2 (pdf)
    output$dlqccorplot2pdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                "dlqccorplot2pdfimg",
                "Download Static R Plot (PDF)"
            )
        }
    })

    ## COR - DOWNLOAD PLOT - corplot2 (pdf)
    output$dlqccorplot2pdfimg <- downloadHandler(
        filename =  function() {
            paste("cor-scatterplot.pdf")
        },
        content = function(file) {
            pdf(file, width = 8, height = 6.5, onefile = FALSE)
            corScatter(
                s.cor = event_data("plotly_click", source = "corplot"),
                cts.tran = corout()[[2]],
                tran = input$transform,
                lab = ddstran()[[2]]
            )
            dev.off()
        }
    )

    ## COR - DOWNLOAD BUTTON - corplot2 (png)
    output$dlqccorplot2png <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                "dlqccorplot2pngimg",
                "Download Static R Plot (PNG)"
            )
        }
    })

    ## COR - DOWNLOAD PLOT - corplot2 (png)
    output$dlqccorplot2pngimg <- downloadHandler(
        filename =  function() {
            paste("cor-scatterplot.png")
        },
        content = function(file) {
            corScatter(
                s.cor = event_data("plotly_click", source = "corplot"),
                cts.tran = corout()[[2]],
                tran = input$transform,
                lab = ddstran()[[2]]
            )
            dev.off()
        }
    )

    ## COR - header (2) - Correlation analysis
    output$headcor2 <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            h4("Sample Distance Matrix")
        }
    })

    ## COR - visualization - sample distance matrix
    output$corplot3 <- renderPlot({
        if (input$goqc == 0) {
            return()
        } else {
            withProgress(message = "Rendering distance matrix...", value = 0, {
                sampdistPlot(
                    cts = ddstran()[[1]]
                )
            })
        }
    })

    ## COR - Show download button - sample distance matrix (PDF)
    output$dlqcorplot3pdf <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcorplot3pdfimg", "Download Plot (PDF)")
        }
    })

    ## COR - Download plot - sample distance matrix (PDF)
    output$dlqcorplot3pdfimg <- downloadHandler(
        filename =  function() {
            paste("sample-dists.pdf")
        },
        content = function(file) {
            pdf(file, width = 7, height = 6.5, onefile = FALSE)
            sampdistPlot(
                cts = ddstran()[[1]]
            )
            dev.off()
        }
    )

    ## COR - Show download button - sample distance matrix (PNG)
    output$dlqcorplot3png <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton("dlqcorplot3pngimg", "Download Plot (PNG)")
        }
    })

    ## COR - Download plot - sample distance matrix (PNG)
    output$dlqcorplot3pngimg <- downloadHandler(
        filename =  function() {
            paste("sample-dists.png")
        },
        content = function(file) {
            png(file, width = 800, height = 750)
            sampdistPlot(
                cts = ddstran()[[1]]
            )
            dev.off()
        }
    )





    ###################################################################
    ###################################################################
    ### SECTION 07 - SESSION INFO (SESSION)
    ###################################################################
    ###################################################################

    ## SESSION - Add session info verbatim
    output$sessinfo <- renderPrint({
            sessionInfo()
    })





    ###################################################################
    ###################################################################
    ### SECTION 08 - GEO METADATA GENERATION (GEO)
    ###################################################################
    ###################################################################

    #------------------------------------------------------------------
    # GEO - Series
    #------------------------------------------------------------------

    output$geo_series <- renderUI({
        if(input$goqc == 0) {
            p(
                em("Please load data first!"),
                style = "color:grey"
            )
        } else {
            list(
                h4("Series (1 - 6)"),
                p(
                    HTML(
                        paste0(
                            "This section describes the overall experiment. ",
                            "Once finished, click on the ",
                            "<b>Submit Series Info</b> button."
                        )
                    )
                )
            )
        }
    })

    output$geo_01_title <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            textInput(
                inputId = "geo_01_title",
                label = "1. Title",
                value = "",
                width = "500px"
            )
        }
    })

    output$geo_02_summary <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            textAreaInput2(
                inputId = "geo_02_summary",
                label = "2. Summary",
                value = "",
                width = "500px",
                rows = 5,
                resize = "vertical"
            )
        }
    })

    output$geo_xx_overall_design <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            textAreaInput2(
                inputId = "geo_xx_overall_design",
                label = "3. Overall Design",
                value = "",
                width = "500px",
                rows = 5,
                resize = "vertical"
            )
        }
    })

    output$geo_03a_contrib_action <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton("addcontrib", "Add Contributor"),
                actionButton("rmvcontrib", "Remove Contributor")
            )
        }
    })

    values <- reactiveValues(num_contrib = 0)
    observeEvent(input$addcontrib, ignoreNULL = FALSE, {
        if(input$goqc == 0) {
            return()
        } else {
            values$num_contrib <- values$num_contrib + 1
            num <- values$num_contrib
            insertUI(
                selector = "#contrib",
                where = "beforeEnd",
                div(
                    id = paste0("contrib", num),
                    textInput(
                        inputId = paste0(
                            "geo_contrib_",
                            str_pad(num, 3, pad = "0")
                        ),
                        label = paste0(
                            "4", LETTERS[num], ". Contributor ", num
                        ),
                        width = "500px"
                    )
                )
            )
        }
    })



    observeEvent(input$rmvcontrib, {
        num <- values$num_contrib
        # Don't let the user remove the very first contrib
        if (num == 1) {
            return()
        } else {
            removeUI(selector = paste0("#contrib", num))
            values$num_contrib <- values$num_contrib - 1
            updateTextInput(
                session = session,
                inputId = paste0(
                    "geo_contrib_",
                    str_pad(num, 3, pad = "0")
                ),
                label = paste0(
                    "4", LETTERS[num], ". Contributor ", num
                ),
                value = "NA"
            )
        }
    })


    output$geo_04_suppfile <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                br(),
                br(),
                textInput(
                    inputId = "geo_04_suppfile",
                    label = "5. Supplementary File (optional)",
                    value = "",
                    width = "500px"
                )
            )
        }
    })

    output$geo_05_sra <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                textInput(
                    inputId = "geo_05_sra",
                    label = "6. SRA Center Name Code (optional)",
                    value = "",
                    width = "500px"
                )
            )
        }
    })

    output$geo_submit_series <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_series",
                    label = "Submit Series Info",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_series_check <- renderUI({
        validate(
            need(input$geo_submit_series != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_01 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                hr(),
                br()
            )
        }
    })


    #------------------------------------------------------------------
    # GEO - Samples
    #------------------------------------------------------------------

    output$geo_samples <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                h4("Samples (7)"),
                p(
                    paste(
                        "This section lists and describes each of the",
                        "biological Samples under investgation, as well as",
                        "any protocols that are specific to individual",
                        "Samples. Additional \"processed data files\" or",
                        "\"raw files\" may be included."
                    )
                ),
                p(
                    HTML(
                        paste0(
                            "Based on your metadata, you have ",
                            "<b>", nrow(ddsout()[[2]]), "</b> ",
                            "samples. Please fill out the ",
                            "following tabs. Once finished, click on the ",
                            "<b>Submit Sample Info</b> button."
                        )

                    )
                )
            )
        }
    })

    output$geo_samples_list <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            nTabs <- nrow(ddsout()[[2]])

            myTabs <- lapply(seq_len(nTabs), function(i) {
                tabPanel(
                    title = paste0("Sample ", i),
                    br(),
                    h5(
                        HTML(
                            paste0(
                                "Based on your metadata, Sample ", i, " is ",
                                "\"",
                                "<b>", rownames(ddsout()[[2]])[i], "</b>",
                                "\".",
                                " You may change the title section if you",
                                " want something more descriptive."
                            )
                        )
                    ),
                    br(),
                    uiOutput(paste0("geo_sample_title_", i)),
                    uiOutput(paste0("geo_sample_source_", i)),
                    uiOutput(paste0("geo_sample_organism_", i)),
                    uiOutput(paste0("geo_sample_molecule_", i)),
                    uiOutput(paste0("geo_sample_description_", i)),
                    br(),
                    p(
                        HTML(
                            paste0(
                                "<b> 7G", i, "a. Processed Data File 1</b>"
                            )
                        )
                    ),
                    p(
                        HTML(
                            paste0(
                                "<b>Note:</b> ",
                                "Your processed data file name for the ",
                                "sample ", "\"",
                                "<b>", rownames(ddsout()[[2]])[i], "</b>",
                                "\"", " is ", "\"",
                                "<b>", rownames(ddsout()[[2]])[i],
                                ".txt </b>", "\".",
                                "This file is the raw counts pertaining ",
                                "to that sample. ",
                                "If you would like to ",
                                "add more processed data files, please enter ",
                                "values below. If not, <i>please leave blank",
                                "</i>!"
                            )
                        )
                    ),
                    div(id = paste0("geo_pdfile_", i)),
                    uiOutput(paste0("geo_pdfile_action_", i)),
                    br(),
                    br(),
                    div(id = paste0("geo_rawfile_", i)),
                    uiOutput(paste0("geo_rawfile_action_", i)),
                    br()
                )
            })
            do.call(tabsetPanel, myTabs)
        }
    })

    output$geo_submit_sample <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_sample",
                    label = "Submit Sample Info",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_sample_check <- renderUI({
        validate(
            need(input$geo_submit_sample != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_02 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                hr(),
                br()
            )
        }
    })

    observe(
        lapply(seq_len(nrow(ddsout()[[2]])), function(i) {
            output[[paste0("geo_sample_title_", i)]] <- renderUI({
                textInput(
                    inputId = paste0(
                        "geo_sample_title_",
                        str_pad(i, 3, pad = "0")
                    ),
                    label = paste0("7A", i, ".", " Title"),
                    value = paste0(
                        rownames(ddsout()[[2]])[i]
                    ),
                    width = "500px"
                )
            })

            output[[paste0("geo_sample_source_", i)]] <- renderUI({
                textInput(
                    inputId = paste0(
                        "geo_sample_source_",
                        str_pad(i, 3, pad = "0")
                    ),
                    label = paste0("7B", i, ".", " Source Name"),
                    value = "",
                    width = "500px"
                )
            })

            output[[paste0("geo_sample_organism_", i)]] <- renderUI({
                textInput(
                    inputId = paste0(
                        "geo_sample_organism_",
                        str_pad(i, 3, pad = "0")
                    ),
                    label = paste0("7C", i, ".", " Organism"),
                    value = "",
                    width = "500px"
                )
            })

            output[[paste0("geo_sample_molecule_", i)]] <- renderUI({
                selectInput(
                    inputId = paste0("geo_sample_molecule_", i),
                    label = paste0("7D", i, ".", " Molecule"),
                    choices = c(
                        "total RNA" = "total RNA",
                        "polyA RNA" = "polyA RNA",
                        "cytoplasmic RNA" = "cytoplasmic RNA",
                        "nuclear RNA" = "nuclear RNA",
                        "genomic DNA" = "genomic DNA",
                        "protein" = "protein",
                        "other" = "other"
                    ),
                    width = "500px"
                )
            })

            output[[paste0("geo_sample_description_", i)]] <- renderUI({
                textInput(
                    inputId = paste0("geo_sample_description_", i),
                    label = paste0("7E", i, ".", " Description"),
                    value = "",
                    width = "500px"
                )
            })

            output[[paste0("geo_pdfile_action_", i)]] <- renderUI({
                list(
                    actionButton(
                        paste0("add_geo_pdfile_", i),
                        "Add Proc. Data File"
                    ),
                    actionButton(
                        paste0("rmv_geo_pdfile_", i),
                        "Remove Proc. Data File"
                    )
                )
            })

            values3 <- reactiveValues(num_pdfile = -1)
            observeEvent(
                input[[paste0("add_geo_pdfile_", i)]],
                ignoreNULL = FALSE, {
                if (input$goqc == 0) {
                    return()
                } else {
                    values3$num_pdfile <- values3$num_pdfile + 1
                    num <- values3$num_pdfile
                    insertUI(
                        selector = paste0("#geo_pdfile_", i),
                        where = "beforeEnd",
                        div(
                            id = paste0("geo_pdfile_", i, num),
                            textInput(
                                inputId = paste0(
                                    "geo_pdfile_",
                                    str_pad(i, 3, pad = "0"),
                                    "_",
                                    str_pad(num + 1, 3, pad = "0")
                                ),
                                label = paste0(
                                    "7F", i, letters[num + 1],
                                    ". Processed Data File ", num + 1
                                ),
                                width = "500px"
                            )
                        )
                    )
                }
            })

            observeEvent(input[[paste0("rmv_geo_pdfile_", i)]], {
                num <- values3$num_pdfile
                # Don't let the user remove the very first contrib
                if (num == 1) {
                    return()
                } else {
                    removeUI(
                        selector = paste0("#geo_pdfile_", i, num)
                    )
                    values3$num_pdfile <- values3$num_pdfile - 1
                    updateTextInput(
                        session = session,
                        inputId = paste0(
                            "geo_pdfile_",
                            str_pad(i, 3, pad = "0"),
                            "_",
                            str_pad(num + 1, 3, pad = "0")
                        ),
                        label = paste0(
                            "7F", i, letters[num + 1],
                            ". Processed Data File ", num + 1
                        ),
                        value = "NA"
                    )
                }
            })

            output[[paste0("geo_rawfile_action_", i)]] <- renderUI({
                list(
                    actionButton(
                        paste0("add_geo_rawfile_", i),
                        "Add Raw Data File"
                    ),
                    actionButton(
                        paste0("rmv_geo_rawfile_", i),
                        "Remove Raw Data File"
                    )
                )
            })

            values4 <- reactiveValues(num_rawfile = -1)
            observeEvent(
                input[[paste0("add_geo_rawfile_", i)]],
                ignoreNULL = FALSE, {
                if (input$goqc == 0) {
                    return()
                } else {
                    values4$num_rawfile <- values4$num_rawfile + 1
                    num <- values4$num_rawfile
                    insertUI(
                        selector = paste0("#geo_rawfile_", i),
                        where = "beforeEnd",
                        div(
                            id = paste0("geo_rawfile_", i, num),
                            textInput(
                                inputId = paste0(
                                    "geo_rawfile_",
                                    str_pad(i, 3, pad = "0"),
                                    "_",
                                    str_pad(num, 3, pad = "0")
                                ),
                                label = paste0(
                                    "7G", i, letters[num],
                                    ". Raw Data File ", num
                                ),
                                width = "500px"
                            )
                        )
                    )
                }
            })

            observeEvent(input[[paste0("rmv_geo_rawfile_", i)]], {
                num <- values4$num_rawfile
                # Don't let the user remove the very first contrib
                if (num == 1) {
                    return()
                } else {
                    removeUI(selector = paste0("#geo_rawfile_", i, num))
                    values4$num_rawfile <- values4$num_rawfile - 1
                    updateTextInput(
                        session = session,
                        inputId = paste0(
                            "geo_rawfile_",
                            str_pad(i, 3, pad = "0"),
                            "_",
                            str_pad(num, 3, pad = "0")
                        ),
                        label = paste0(
                            "7G", i, letters[num],
                            ". Raw Data File ", num
                        ),
                        value = ""
                    )
                }
            })
        })
    )



    #------------------------------------------------------------------
    # GEO - Protocols
    #------------------------------------------------------------------

    output$geo_protocols <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                h4("Protocols (8 - 12)"),
                p(
                    HTML(
                        paste(
                            "Please provide information about what",
                            "protocols were used to conduct this experiment."
                        )

                    )
                )
            )
        }
    })

    output$geo_07_growth_protocol <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            textAreaInput2(
                inputId = "geo_07_growth_protocol",
                label = "8. Growth Protocol (optional)",
                value = "",
                width = "500px",
                rows = 5,
                resize = "vertical"
            )
        }
    })

    output$geo_08_treatment_protocol <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                textAreaInput2(
                    inputId = "geo_08_treatment_protocol",
                    label = "9. Treatment Protocol (optional)",
                    value = "",
                    width = "500px",
                    rows = 5,
                    resize = "vertical"
                )
            )
        }
    })

    output$geo_09_extract_protocol <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                textAreaInput2(
                    inputId = "geo_09_extract_protocol",
                    label = "10. Extract Protocol",
                    value = "",
                    width = "500px",
                    rows = 5,
                    resize = "vertical"
                )
            )
        }
    })

    output$geo_10_lib_construct_protocol <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                textAreaInput2(
                    inputId = "geo_10_lib_construct_protocol",
                    label = "11. Libary Construction Protocol",
                    value = "",
                    width = "500px",
                    rows = 5,
                    resize = "vertical"
                )
            )
        }
    })

    output$geo_11_lib_strategy <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                selectInput(
                    inputId = "geo_11_lib_strategy",
                    label = "12. Library Strategy",
                    choices = c(
                        "RNA-Seq" = "RNA-Seq",
                        "miRNA-Seq" = "miRNA-Seq",
                        "ncRNA-Seq" = "ncRNA-Seq",
                        "RNA-Seq (CAGE)" = "RNA-Seq (CAGE)",
                        "RNA-Seq (RACE)" = "RNA-Seq (RACE)",
                        "ChIP-Seq" = "ChIP-Seq",
                        "MNase-Seq" = "MNase-Seq",
                        "MBD-Seq" = "MBD-Seq",
                        "MRE-Seq" = "MRE-Seq",
                        "Bisulfite-Seq" = "Bisulfite-Seq",
                        "Bisulfite-Seq (reduced representation)" =
                            "Bisulfite-Seq (reduced representation)",
                        "MeDIP-Seq" = "MeDIP-Seq",
                        "DNAse-Hypersensitivity" = "DNAse-Hypersensitivity",
                        "Tn-Seq" = "Tn-Seq",
                        "FAIRE-seq" = "FAIRE-seq",
                        "SELEX" = "SELEX",
                        "RIP-Seq" = "RIP-Seq",
                        "ChIA-PET" = "ChIA-PET",
                        "OTHER" = "OTHER"
                    ),
                    width = "500px"
                )
            )
        }
    })

    output$geo_submit_protocol <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_protocol",
                    label = "Submit Protocol Info",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_protocol_check <- renderUI({
        validate(
            need(input$geo_submit_protocol != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_03 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                hr(),
                br()
            )
        }
    })



    #------------------------------------------------------------------
    # GEO - Data processing pipeline
    #------------------------------------------------------------------

    output$geo_data_proc_pipeline <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                h4("Data Processing Pipeline (13 - 15)"),
                p(
                    HTML(
                        paste(
                            "Data processing steps include base-calling, ",
                            "alignment, filtering, peak-calling, generation",
                            "of normalized abundance measurements, etc. For",
                            "each step, provide a description, as well as",
                            "software name, version, parameters, if",
                            "applicable. Once finished, click on the ",
                            "<b>Submit Pipeline Info</b> button."
                        )

                    )
                )
            )
        }
    })

    output$geo_12a_data_proc_action <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton("adddatastep", "Add Data Proc. Step"),
                actionButton("rmvdatastep", "Remove Data Proc. Step")
            )
        }
    })

    values5 <- reactiveValues(num_data_step = 0)
    observeEvent(input$adddatastep, ignoreNULL = FALSE, {
        if(input$goqc == 0) {
            return()
        } else {
            values5$num_data_step <- values5$num_data_step + 1
            num <- values5$num_data_step
            insertUI(
                selector = "#data_step",
                where = "beforeEnd",
                div(
                    id = paste0("data_step_", num),
                    textInput(
                        inputId = paste0(
                            "data_step_",
                            str_pad(num, 3, pad = "0")
                        ),
                        label = paste0(
                            "13", LETTERS[num],
                            ". Data Processing Step ", num
                        ),
                        width = "500px"
                    )
                )
            )
        }
    })

    observeEvent(input$rmvdatastep, {
        num <- values5$num_data_step
        if (num == 1) {
            return()
        } else {
            removeUI(selector = paste0("#data_step_", num))
            values5$num_data_step <- values5$num_data_step - 1
            updateTextInput(
                session = session,
                inputId = paste0(
                    "data_step_",
                    str_pad(num, 3, pad = "0")
                ),
                label = paste0(
                    "13", LETTERS[num],
                    ". Data Processing Step ", num
                ),
                value = "NA"
            )
        }
        # removeUI(selector = paste0("#data_step_", num))
        # values5$num_data_step <- values5$num_data_step - 1
    })

    output$geo_13_genome_build <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                br(),
                br(),
                textInput(
                    inputId = "geo_13_genome_build",
                    label = "14. Genome Build",
                    value = "",
                    width = "500px"
                )
            )
        }
    })

    output$geo_14_proc_data_files <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                textInput(
                    inputId = "geo_14_proc_data_files",
                    label = "15. Processed Data Files Format and Content",
                    value = "",
                    width = "500px"
                )
            )
        }
    })

    output$geo_submit_pipeline <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_pipeline",
                    label = "Submit Pipeline Info",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_pipeline_check <- renderUI({
        validate(
            need(input$geo_submit_pipeline != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_04 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                hr(),
                br()
            )
        }
    })



    #------------------------------------------------------------------
    # GEO - processed data
    #------------------------------------------------------------------

    output$geo_proc_head <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            list(
                h4("Processed Data Files (16)")
            )
        }
    })


    output$geo_proc_data <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(
                input$geo_submit_sample != 0,
                paste(
                    "You have not submitted your sample info yet.",
                    "This tab section will not be populated until",
                    "you have entered additional processed data",
                    "via the \"Samples\" section."
                )
            )
        )

        x <- pdf_out()
        if (length(x) == 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "Based on your sample information, you have",
                            "<b>", length(x), "</b>",
                            "additional processed data files.",
                            "This will <b>not</b> affect your final metadata",
                            "file."
                        )
                    )
                )
            )
        } else {
            list(
                p(
                    paste(
                        "For each file listed in the \"processed data",
                        "file\" columns of the SAMPLES section, provide",
                        "additional information below."
                    )
                ),
                p(
                    HTML(
                        paste(
                            "Based on your sample information, you have",
                            "<b>", length(x), "</b>",
                            "addtional processed data file(s).",
                            "Please fill out the following tabs.",
                            "Once finished, click on the ",
                            "<b>Submit Additional Processed Data</b> button."
                        )
                    )
                )
            )
        }
    })

    output$geo_proc_data_list <- renderUI({
        x <- pdf_out()

        if(length(x) == 0) {
            return()
        } else {
            myTabs <- lapply(seq_len(length(x)), function(i) {
                tabPanel(
                    paste0("Proc. Data File ", i),
                    br(),
                    h5(
                        HTML(
                            paste0(
                                "Based on your sample submission, Proc. Data",
                                " File ", i, " is ",
                                "\"", "<b>", paste0(x[[i]]), "</b>", "\"."
                            )
                        )
                    ),
                    br(),
                    textInput(
                        inputId = paste0("filetype_", names(x)[i]),
                        label = paste0("16A", i, ". File Type"),
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0("filecheck_", names(x)[i]),
                        label = paste0("16B", i, ". MD5 File Checksum"),
                        width = "500px"
                    )
                )
            })
            do.call(tabsetPanel, myTabs)
        }
    })

    output$geo_submit_proc_data <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(input$geo_submit_sample != 0, "")
        )
        x <- pdf_out()

        if (length(x) == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_proc_data",
                    label = "Submit Additional Processed Data",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_proc_data_check <- renderUI({
        validate(
            need(input$geo_submit_proc_data != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_05 <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        list(
            hr(),
            br()
        )
    })



    #------------------------------------------------------------------
    # GEO - Raw data
    #------------------------------------------------------------------

    output$geo_raw_head <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            list(
                h4("Raw Data Files (17)")
            )
        }
    })

    output$geo_raw_data <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(
                input$geo_submit_sample != 0,
                paste(
                    "You have not submitted your sample info yet.",
                    "This tab section will not be populated until",
                    "you have entered raw data information."
                )
            )
        )
        x <- raw_out()

        if (length(x) == 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "You have not yet entered sample information",
                            "pertaining to raw data files.",
                            "<b>Please enter this info in order",
                            "to continue!</b>"
                        )
                    )
                )
            )
        } else {
            list(
                p(
                    HTML(
                        paste(
                                "For each file listed in the \"raw",
                                "file\" columns of the SAMPLES section,",
                                "provide additional information below.",
                                "Once finished, click on the ",
                                "<b>Submit Raw Data</b> button."
                        )

                    )
                ),
                p(
                    HTML(
                        paste(
                            "Based on your sample information, you have",
                            "<b>", length(x), "</b>",
                            "raw data file(s).",
                            "Please fill out the following tabs."
                        )

                    )
                )
            )
        }
    })

    output$geo_raw_data_list <- renderUI({
        x <- raw_out()
        if (length(x) == 0) {
            return()
        } else {
            myTabs <- lapply(seq_len(length(x)), function(i) {
                tabPanel(
                    paste0("Raw Data File ", i),
                    br(),
                    h5(
                        HTML(
                            paste0(
                                "Based on your submission, Raw Data File ",
                                i, " is ",
                                "\"", "<b>",
                                paste0(x[[i]]), "</b>", "\". "
                            )
                        )
                    ),
                    br(),
                    textInput(
                        inputId = paste0("filetype_", names(x)[i]),
                        label = paste0("17A", i, ". File Type"),
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0("filecheck_", names(x)[i]),
                        label = paste0("17B", i, ". MD5 File Checksum"),
                        width = "500px"
                    ),
                    selectInput(
                        inputId = paste0("instmod_", names(x)[i]),
                        label = paste0("17C", i, ". Instrument Model"),
                        choices = c(
                            "Illumina Genome Analyzer" =
                                "Illumina Genome Analyzer",
                            "Illumina Genome Analyzer II" =
                                "Illumina Genome Analyzer II",
                            "Illumina Genome Analyzer IIx" =
                                "Illumina Genome Analyzer IIx",
                            "Illumina HiSeq 2000" =
                                "Illumina HiSeq 2000",
                            "Illumina HiSeq 1000" =
                                "Illumina HiSeq 1000",
                            "Illumina MiSeq" =
                                "Illumina MiSeq",
                            "AB SOLiD System" =
                                "AB SOLiD System",
                            "AB SOLiD System 2.0" =
                                "AB SOLiD System 2.0",
                            "AB SOLiD System 3.0" =
                                "AB SOLiD System 3.0",
                            "AB SOLiD 4 System" =
                                "AB SOLiD 4 System",
                            "AB SOLiD 4hq System" =
                                "AB SOLiD 4hq System",
                            "AB SOLiD PI System" =
                                "AB SOLiD PI System",
                            "AB 5500xl Genetic Analyzer" =
                                "AB 5500xl Genetic Analyzer",
                            "AB 5500 Genetic Analyzer" =
                                "AB 5500 Genetic Analyzer",
                            "454 GS" = "454 GS",
                            "454 GS 20" = "454 GS 20",
                            "454 GS FLX" = "454 GS FLX",
                            "454 GS Junior" = "454 GS Junior",
                            "454 GS FLX Titanium" = "454 GS FLX Titanium",
                            "Helicos HeliScope" = "Helicos HeliScope",
                            "PacBio RS" = "PacBio RS",
                            "PacBio RS II" = "PacBio RS II",
                            "Complete Genomics" = "Complete Genomics",
                            "Ion Torrent PGM" = "Ion Torrent PGM"
                        ),
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0("readlen_", names(x)[i]),
                        label = paste0("17D", i, ". Read Length"),
                        width = "500px"
                    ),
                    radioButtons(
                        inputId = paste0("spend_", names(x)[i]),
                        label = paste0(
                            "16E", i, ". Single, Paired-end or SOLiD?"
                        ),
                        choices = c(
                            "Single" = "single",
                            "Paired-end" = "paired-end",
                            "SOLiD" = "paired-end (SOLiD)"
                        )
                    )
                )
            })
        }
        do.call(tabsetPanel, myTabs)
    })

    output$geo_submit_raw_data <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(input$geo_submit_sample != 0, "")
        )
        x <- raw_out()

        if (length(x) == 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_raw_data",
                    label = "Submit Raw Data",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_raw_data_check <- renderUI({
        validate(
            need(input$geo_submit_raw_data != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_06 <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        list(
            hr(),
            br()
        )
    })



    #------------------------------------------------------------------
    # GEO - paired end data
    #------------------------------------------------------------------

    output$geo_pair_head <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            list(
                h4("Paired-end Experiments (18)")
            )
        }
    })

    output$geo_paired_end <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(
                input$geo_submit_raw_data != 0,
                paste(
                    "You have not submitted your raw data info yet.",
                    "This tab section will not be populated until",
                    "you have submitted single read/paired-end info in the",
                    "\"Raw Data Files\" section."
                )
            )
        )
        x <- raw_data_file_out()
        x_pair <- x[grep("^spend_geo_rawfile_", names(x))]
        x_pair <- x_pair[x_pair == "paired-end"]

        if (length(x_pair) == 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "Your raw data does not contain any paired-end",
                            "data. This will <b>not</b> affect your",
                            "final metadata file."
                        )
                    )
                )
            )
        } else if (length(x_pair) %% 2 != 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "It seems that you have an uneven number of",
                            "paired-end read files. Please make sure you have",
                            "<i>two reads per sample.</i>"
                        )
                    )
                )
            )
        } else {
            list(
                p(
                    HTML(
                        paste(
                            "For paired-end experiments, provide the",
                            "average insert size and standard deviation",
                            "if known."
                        )

                    )
                ),
                p(
                    HTML(
                        paste(
                            "Based on your sample information, you have",
                            "<b>", length(x_pair), "</b>",
                            "file(s) that are paired-end.",
                            "Please fill out the following tab(s)."
                        )

                    )
                )
            )
        }
    })

    output$geo_paired_end_list <- renderUI({
        # Get raw data metadata...
        x <- raw_data_file_out()
        x_pair <- x[grep("^spend_geo_rawfile_", names(x))]
        x_pair <- x_pair[x_pair == "paired-end"]

        # Get all raw data
        y <- raw_out()

        # Get raw data that's paired-end
        tmp_x <- x_pair
        names(tmp_x) <- gsub("^spend_", "", names(tmp_x))
        tmp_y <- y[intersect(names(tmp_x), names(y))]
        tmp_y <- unlist(tmp_y, use.names = FALSE)

        x_inst <- x[grep("^instmod_geo_rawfile_", names(x))]
        x_inst_solid <- x_inst[x_inst == "AB SOLiD System"
        | x_inst == "AB SOLiD System 2.0"
        | x_inst == "AB SOLiD System 3.0"
        | x_inst == "AB SOLiD 4 System"
        | x_inst == "AB SOLiD 4hq System"
        | x_inst == "AB SOLiD PI System"]

        if(length(x_pair) == 0 | length(x_pair) %% 2 != 0) {
            return()
        } else if (length(x_inst_solid) == 0) {
            myTabs <- lapply(seq_len(length(x_pair) / 2), function(i) {
                tabPanel(
                    paste0("Paired-end Read Set ", i),
                    br(),
                    h5(
                        HTML(
                            paste0(
                                "Please enter which files are <b>Read 1</b> ",
                                "and <b>Read 2</b>:"
                            )
                        )
                    ),
                    br(),
                    selectInput(
                        inputId = paste0("r1_row_", str_pad(i, 3, pad = "0")),
                        label = paste0("18A", i, ". File Name 1"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    selectInput(
                        inputId = paste0("r2_row_", str_pad(i, 3, pad = "0")),
                        label = paste0("18B", i, ". File Name 2"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0("ais_", str_pad(i, 3, pad = "0")),
                        label = paste0("18C", i, ". Average Insert Size"),
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0("std_", str_pad(i, 3, pad = "0")),
                        label = paste0("18D", i, ". Standard Deviation"),
                        width = "500px"
                    )
                )
            })
            do.call(tabsetPanel, myTabs)
        }
    })


    output$geo_submit_paired_end <- renderUI({
        # Get raw data metadata...
        x <- raw_data_file_out()
        x_pair <- x[grep("^spend_geo_rawfile_", names(x))]
        x_pair <- x_pair[x_pair == "paired-end"]

        # Get all raw data
        y <- raw_out()

        # Get raw data that's paired-end
        tmp_x <- x_pair
        names(tmp_x) <- gsub("^spend_", "", names(tmp_x))
        tmp_y <- y[intersect(names(tmp_x), names(y))]
        tmp_y <- unlist(tmp_y, use.names = FALSE)

        x_inst <- x[grep("^instmod_geo_rawfile_", names(x))]
        x_inst_solid <- x_inst[x_inst == "AB SOLiD System"
        | x_inst == "AB SOLiD System 2.0"
        | x_inst == "AB SOLiD System 3.0"
        | x_inst == "AB SOLiD 4 System"
        | x_inst == "AB SOLiD 4hq System"
        | x_inst == "AB SOLiD PI System"]

        if(length(x_pair) == 0 | length(x_pair) %% 2 != 0) {
            return()
        } else if (length(x_inst_solid) == 0) {
            list(
                actionButton(
                    inputId = "geo_submit_paired_end",
                    label = "Submit Paired-end Data",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_paired_end_check <- renderUI({
        validate(
            need(input$geo_submit_paired_end != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_07 <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        list(
            hr(),
            br()
        )
    })



    #------------------------------------------------------------------
    # GEO - paired end data - SOLiD
    #------------------------------------------------------------------

    output$geo_solid_head <- renderUI({
        if (input$goqc == 0) {
            return()
        } else {
            list(
                h4("Paired-end Experiments - SOLiD Data (19)")
            )
        }
    })

    output$geo_solid <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        validate(
            need(
                input$geo_submit_raw_data != 0,
                paste(
                    "You have not submitted your raw data info yet.",
                    "This tab section will not be populated until",
                    "you have submitted single read/paired-end info in the",
                    "\"Raw Data Files\" section."
                )
            )
        )
        x <- raw_data_file_out()
        x_solid <- x[grep("^spend_geo_rawfile_", names(x))]
        x_solid <- x_solid[x_solid == "paired-end (SOLiD)"]

        if (length(x_solid) == 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "Your raw data does not contain any SOLiD",
                            "experimental data. This will <b>not</b>",
                            "affect your final metadata file."
                        )
                    )
                )
            )
        } else if (length(x_solid) %% 4 != 0) {
            list(
                p(
                    HTML(
                        paste(
                            "<b>Note:</b>",
                            "It seems that you have an uneven number of",
                            "SOLiD read files. Please make sure you have",
                            "<i>four reads per sample.</i>"
                        )
                    )
                )
            )
        } else {
            list(
                p(
                    HTML(
                        paste(
                            "For SOLiD paired-end experiments, provide the",
                            "average insert size and standard deviation",
                            "if known."
                        )

                    )
                ),
                p(
                    HTML(
                        paste(
                            "Based on your sample information, you have",
                            "<b>", length(x_solid), "</b>",
                            "files that use SOLiD technology.",
                            "Please fill out the following tab(s)."
                        )

                    )
                )
            )
        }
    })

    output$geo_solid_list <- renderUI({
        # Get raw data metadata...
        x <- raw_data_file_out()
        x_solid <- x[grep("^spend_geo_rawfile_", names(x))]
        x_solid <- x_solid[x_solid == "paired-end (SOLiD)"]

        # Get all raw data
        y <- raw_out()

        # Get raw data that's paired-end
        tmp_x <- x_solid
        names(tmp_x) <- gsub("^spend_", "", names(tmp_x))
        tmp_y <- y[intersect(names(tmp_x), names(y))]
        tmp_y <- unlist(tmp_y, use.names = FALSE)

        # Generate tabs
        if(length(x_solid) == 0 | length(x_solid) %% 4 != 0) {
            return()
        } else {
            myTabs <- lapply(seq_len(length(x_solid) / 4), function(i) {
                tabPanel(
                    paste0("SOLiD Set ", i),
                    br(),
                    h5(
                        HTML(
                            paste0(
                                "Please enter which files are <b>Read 1,</b> ",
                                "<b>Read 2</b>, <b>Read 3</b>, and ",
                                "<b>Read 4</b>:"
                            )
                        )
                    ),
                    br(),
                    selectInput(
                        inputId = paste0(
                            "solid_r1_row_", str_pad(i, 3, pad = "0")
                        ),
                        label = paste0("19A", i, ". File Name 1"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    selectInput(
                        inputId = paste0(
                            "solid_r2_row_", str_pad(i, 3, pad = "0")
                        ),
                        label = paste0("19B", i, ". File Name 2"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    selectInput(
                        inputId = paste0(
                            "solid_r3_row_", str_pad(i, 3, pad = "0")
                        ),
                        label = paste0("19C", i, ". File Name 3"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    selectInput(
                        inputId = paste0(
                            "solid_r4_row_", str_pad(i, 3, pad = "0")
                        ),
                        label = paste0("19D", i, ". File Name 4"),
                        choices = tmp_y,
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0(
                            "solid_ais_", str_pad(i, 3, pad = "0")),
                        label = paste0("19E", i, ". Average Insert Size"),
                        width = "500px"
                    ),
                    textInput(
                        inputId = paste0(
                            "solid_std_", str_pad(i, 3, pad = "0")
                        ),
                        label = paste0("19F", i, ". Standard Deviation"),
                        width = "500px"
                    )
                )
            })
            do.call(tabsetPanel, myTabs)
        }
    })

    output$geo_submit_solid <- renderUI({
        # Get raw data metadata...
        x <- raw_data_file_out()
        x_solid <- x[grep("^spend_geo_rawfile_", names(x))]
        x_solid <- x_solid[x_solid == "paired-end (SOLiD)"]

        # Get all raw data
        y <- raw_out()

        # Get raw data that's paired-end
        tmp_x <- x_solid
        names(tmp_x) <- gsub("^spend_", "", names(tmp_x))
        tmp_y <- y[intersect(names(tmp_x), names(y))]
        tmp_y <- unlist(tmp_y, use.names = FALSE)

        if(length(x_solid) == 0 | length(x_solid) %% 4 != 0) {
            return()
        } else {
            list(
                actionButton(
                    inputId = "geo_submit_solid",
                    label = "Submit SOLiD Data",
                    icon = icon("table")
                )
            )
        }
    })

    output$geo_submit_solid_check <- renderUI({
        validate(
            need(input$geo_submit_solid != 0, "")
        )
        list(
            h5(
                HTML("<b>Submitted!</b>")
            )
        )
    })

    output$geo_space_08 <- renderUI({
        validate(
            need(input$goqc != 0, "")
        )
        list(
            hr(),
            br()
        )
    })



    #------------------------------------------------------------------
    # GEO - Download options
    #------------------------------------------------------------------

    geo_ser_out <- eventReactive(input$geo_submit_series, {
        x <- reactiveValuesToList(input)
        x1 <- x[grep("geo_01_title", names(x))]
        x2 <- x[grep("geo_02_summary", names(x))]
        xX <- x[grep("geo_xx_overall_design", names(x))]
        x3 <- x[grep("^geo_contrib_", names(x))]
        x3 <- x3[!x3 %in% "NA"]
        x4 <- x[grep("geo_04_suppfile", names(x))]
        x5 <- x[grep("geo_05_sra", names(x))]
        x <- c(x1, x2, xX, x3, x4, x5)
        return(x)
    })

    geo_sam_out_01 <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x1 <- x[grep("^geo_sample_title_", names(x))]
        x2 <- x[grep("^geo_sample_source_", names(x))]
        x3 <- x[grep("^geo_sample_organism_", names(x))]
        return(list(x1, x2, x3))
    })

    # geo_sam_out_02 <- eventReactive(input$geo_submit_sample, {
    #     x <- reactiveValuesToList(input)
    #     x <- x[grep("^geo_character_", names(x))]
    #     x <- x[!x %in% "NA"]
    #     x <- x[order(names(x))]
    #     x <- char_sorter(x, rownames(ddsout()[[2]]))
    #     return(x)
    # })

    geo_sam_out_03 <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x1 <- x[grep("^geo_sample_molecule_", names(x))]
        x2 <- x[grep("^geo_sample_description_", names(x))]
        return(list(x1, x2))
    })

    geo_sam_out_04 <- eventReactive(input$geo_submit_sample, {
        x <- rownames(ddsout()[[2]])
        return(x)
    })

    geo_sam_out_05 <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_pdfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[order(names(x))]
        x <- pdf_sorter(x, rownames(ddsout()[[2]]))
        return(x)
    })

    geo_sam_out_06 <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[order(names(x))]
        x <- raw_sorter(x, rownames(ddsout()[[2]]))
        return(x)
    })

    geo_prot_out_01 <- eventReactive(input$geo_submit_protocol, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_0[7-9]|1[0-1]_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_pipe_out_01 <- eventReactive(input$geo_submit_pipeline, {
        x <- reactiveValuesToList(input)
        x1 <- x[grep("^data_step_", names(x))]
        x1 <- x1[!x1 %in% "NA"]
        x2 <- x[grep("^geo_13_genome_build", names(x))]
        x3 <- x[grep("^geo_14_proc_data_files", names(x))]
        x <- c(x1, x2, x3)
        return(x)
    })

    geo_proc_out_01 <- eventReactive(input$geo_submit_proc_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_pdfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_proc_out_02 <- eventReactive(input$geo_submit_proc_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^filetype_geo_pdfile", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_proc_out_03 <- eventReactive(input$geo_submit_proc_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^filecheck_geo_pdfile", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_01 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_02 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^filetype_geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_03 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^filecheck_geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_04 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^instmod_geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_05 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^readlen_geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    geo_raw_out_06 <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^spend_geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    md5_out <- eventReactive(input$geo_submit_sample, {
        x <- cts_out()[[1]]
        return(x)
    })

    cts_out <- eventReactive(input$geo_submit_sample, {
        cts <- ddsout()[[4]]
        cts <- cts_split(cts = cts)
        md5 <- cts[[1]]
        wd <- cts[[2]]
        return(list(md5, wd))
    })

    pdf_out <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_pdfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    raw_out <- eventReactive(input$geo_submit_sample, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^geo_rawfile_", names(x))]
        x <- x[!x %in% "NA"]
        x <- x[!x %in% ""]
        x <- x[order(names(x))]
        return(x)
    })

    raw_data_file_out <- eventReactive(input$geo_submit_raw_data, {
        x <- reactiveValuesToList(input)
        x1 <- x[grep("^filetype_geo_rawfile_", names(x))]
        x2 <- x[grep("^filecheck_geo_rawfile_", names(x))]
        x3 <- x[grep("^instmod_geo_rawfile_", names(x))]
        x4 <- x[grep("^readlen_geo_rawfile_", names(x))]
        x5 <- x[grep("^spend_geo_rawfile_", names(x))]
        x <- c(x1, x2, x3, x4, x5)
        x <- x[order(names(x))]
        return(x)
    })

    geo_pair_out_01 <- eventReactive(input$geo_submit_paired_end, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^r1_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_pair_out_02 <- eventReactive(input$geo_submit_paired_end, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^r2_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_pair_out_03 <- eventReactive(input$geo_submit_paired_end, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^ais_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_pair_out_04 <- eventReactive(input$geo_submit_paired_end, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^std_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_01 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_r1_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_02 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_r2_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_03 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_r3_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_04 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_r4_row_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_05 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_ais_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    geo_solid_out_06 <- eventReactive(input$geo_submit_solid, {
        x <- reactiveValuesToList(input)
        x <- x[grep("^solid_std_", names(x))]
        x <- x[order(names(x))]
        return(x)
    })

    output$geo_debug <- renderPrint({

    })

    output$geo_downloads <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            list(
                h4("Download Options"),
                p(
                    HTML(
                        paste(
                            "Once you have all of your information submitted,",
                            "please click on one (<i>or more</i>) of the",
                            "options: <ul><li><b>Download Metadata</b>",
                            "will return",
                            "the metadata file populated with the prior",
                            "questionnaire in an Excel file format.</li>",
                            "<li><b>Download Processed Data</b> will return a",
                            ".zip file of raw counts for each sample that",
                            "you have submitted to IRIS.</li></ul>",
                            "</br>",
                            "<b>Note</b>: If you are having issues with",
                            "downloading the data, please try the local",
                            "version of the web server. More instructions",
                            "can be found",
                            "<a href=\"https://github.com/btmonier/iris#get-local-application\">here.</a>",
                            "</br></br>"
                        )
                    )
                )
            )
        }
    })

    output$geo_download_excel_1 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                outputId = "geo_download_excel_2",
                label = "Download Metadata (.xlsx)"
            )
        }
    })

    output$geo_download_excel_2 <- downloadHandler(
        filename =  function() {
            "metadata.xlsx"
        },
        content = function(file) {
            wd <- cts_out()[[2]]
            ser_out <- geo_ser_out()
            pdf_out <- pdf_out()
            raw_out <- raw_out()
            sam_title <- rownames(ddsout()[[2]])
            sam_out_01 <- geo_sam_out_01()
            sam_out_02 <- ddsout()[[2]]
            sam_out_03 <- geo_sam_out_03()
            sam_out_04 <- geo_sam_out_04()
            sam_out_05 <- geo_sam_out_05()
            sam_out_06 <- geo_sam_out_06()
            prot_out_01 <- geo_prot_out_01()
            pipe_out_01 <- geo_pipe_out_01()
            md5_out <- md5_out()

            if (length(pdf_out) == 0) {
                proc_out_01 <- list()
                proc_out_02 <- list()
                proc_out_03 <- list()
            } else {
                proc_out_01 <- geo_proc_out_01()
                proc_out_02 <- geo_proc_out_02()
                proc_out_03 <- geo_proc_out_03()
            }

            if (length(raw_out) == 0) {
                raw_out_01 <- list()
                raw_out_02 <- list()
                raw_out_03 <- list()
                raw_out_04 <- list()
                raw_out_05 <- list()
                raw_out_06 <- list()
            } else {
                raw_out_01 <- geo_raw_out_01()
                raw_out_02 <- geo_raw_out_02()
                raw_out_03 <- geo_raw_out_03()
                raw_out_04 <- geo_raw_out_04()
                raw_out_05 <- geo_raw_out_05()
                raw_out_06 <- geo_raw_out_06()
            }

            if (length(raw_out) == 0) {
                pair_out_01 <- list()
                pair_out_02 <- list()
                pair_out_03 <- list()
                pair_out_04 <- list()
                solid_out_01 <- list()
                solid_out_02 <- list()
                solid_out_03 <- list()
                solid_out_04 <- list()
                solid_out_05 <- list()
                solid_out_06 <- list()
            } else if (all(raw_out_06 == "single")) {
                pair_out_01 <- list()
                pair_out_02 <- list()
                pair_out_03 <- list()
                pair_out_04 <- list()
                solid_out_01 <- list()
                solid_out_02 <- list()
                solid_out_03 <- list()
                solid_out_04 <- list()
                solid_out_05 <- list()
                solid_out_06 <- list()
            } else if (all(raw_out_06 == "paired-end (SOLiD)")) {
                pair_out_01 <- list()
                pair_out_02 <- list()
                pair_out_03 <- list()
                pair_out_04 <- list()
                solid_out_01 <- geo_solid_out_01()
                solid_out_02 <- geo_solid_out_02()
                solid_out_03 <- geo_solid_out_03()
                solid_out_04 <- geo_solid_out_04()
                solid_out_05 <- geo_solid_out_05()
                solid_out_06 <- geo_solid_out_06()
            } else if (all(raw_out_06 == "paired-end")) {
                pair_out_01 <- geo_pair_out_01()
                pair_out_02 <- geo_pair_out_02()
                pair_out_03 <- geo_pair_out_03()
                pair_out_04 <- geo_pair_out_04()
                solid_out_01 <- list()
                solid_out_02 <- list()
                solid_out_03 <- list()
                solid_out_04 <- list()
                solid_out_05 <- list()
                solid_out_06 <- list()
            } else if ("paired-end" %in% raw_out_06 & "paired-end (SOLiD)" %in% raw_out_06) {
                pair_out_01 <- geo_pair_out_01()
                pair_out_02 <- geo_pair_out_02()
                pair_out_03 <- geo_pair_out_03()
                pair_out_04 <- geo_pair_out_04()
                solid_out_01 <- geo_solid_out_01()
                solid_out_02 <- geo_solid_out_02()
                solid_out_03 <- geo_solid_out_03()
                solid_out_04 <- geo_solid_out_04()
                solid_out_05 <- geo_solid_out_05()
                solid_out_06 <- geo_solid_out_06()
            } else if ("paired-end" %in% raw_out_06 & !"paired-end (SOLiD)" %in% raw_out_06) {
                pair_out_01 <- geo_pair_out_01()
                pair_out_02 <- geo_pair_out_02()
                pair_out_03 <- geo_pair_out_03()
                pair_out_04 <- geo_pair_out_04()
                solid_out_01 <- list()
                solid_out_02 <- list()
                solid_out_03 <- list()
                solid_out_04 <- list()
                solid_out_05 <- list()
                solid_out_06 <- list()
            } else if (!"paired-end" %in% raw_out_06 & "paired-end (SOLiD)" %in% raw_out_06) {
                pair_out_01 <- list()
                pair_out_02 <- list()
                pair_out_03 <- list()
                pair_out_04 <- list()
                solid_out_01 <- geo_solid_out_01()
                solid_out_02 <- geo_solid_out_02()
                solid_out_03 <- geo_solid_out_03()
                solid_out_04 <- geo_solid_out_04()
                solid_out_05 <- geo_solid_out_05()
                solid_out_06 <- geo_solid_out_06()
            }

            getExcelMetadata(
                wd, ser_out, sam_out_01, sam_out_02, sam_out_03, sam_out_04,
                sam_out_05, sam_out_06, prot_out_01, pipe_out_01, md5_out,
                proc_out_01, proc_out_02, proc_out_03,
                raw_out_01, raw_out_02, raw_out_03,
                raw_out_04, raw_out_05, raw_out_06,
                pair_out_01, pair_out_02, pair_out_03, pair_out_04,
                solid_out_01, solid_out_02, solid_out_03, solid_out_04,
                solid_out_05, solid_out_06
            )

            wd2 <- substring(wd, 5)
            meta <- paste0("metadata-", wd2, ".xlsx")
            file.rename(meta, file)
            setwd("../..")
        },
        contentType = "application/xlsx"
    )

    output$geo_download_proc_data_1 <- renderUI({
        if(input$goqc == 0) {
            return()
        } else {
            downloadButton(
                outputId = "geo_download_proc_data_2",
                label = "Download Processed Data (.zip)"
            )
        }
    })

    output$geo_download_proc_data_2 <- downloadHandler(
        filename = function() {
            "proc-data.zip"
        },
        content = function(file) {
            wd <- cts_out()[[2]]
            setwd(paste0("tmp/", wd, "/"))
            files <- list.files(pattern = "\\.txt$")
            zip(file, files)

            # cts <- ddsout()[[4]]
        }
    )





    ###################################################################
    ###################################################################
    ### SECTION 09 - TAB LINKS (TABS)
    ###################################################################
    ###################################################################

    observeEvent(input$iris_flow_01, {
        updateNavbarPage(
            session,
            "tab_structure",
            selected = "val1"
        )
    }, ignoreInit = TRUE)

    observeEvent(input$iris_flow_02, {
        updateNavbarPage(
            session,
            "tab_structure",
            selected = "val2"
        )
    }, ignoreInit = TRUE)

    observeEvent(input$iris_flow_03, {
        updateNavbarPage(
            session,
            "tab_structure",
            selected = "val3"
        )
    }, ignoreInit = TRUE)

    observeEvent(input$iris_flow_04, {
        updateNavbarPage(
            session,
            "tab_structure",
            selected = "val4"
        )
    }, ignoreInit = TRUE)





    ###################################################################
    ###################################################################
    ### SECTION 10 - BRIC ANALYSIS
    ###################################################################
    ###################################################################

    # DATA - Data load option - submit expression matrix
    output$bric_file1 <- renderUI({
        if (input$bric_examplechoice == "bric_no") {
            fileInput(
                inputId = "bric_file1",
                label = "Submit Expression Matrix",
                accept = c(
                    "text/tab-separated-values"
                )
            )
        } else {
            return()
        }
    })

    # Display cluster number parameter
    output$bric_celltype_number <- renderUI({
        if (input$bric_method == "SC") {
            textInput(
                inputId = "bric_celltype_number",
                label = "Specify number of cell types!"
            )
        } else {
            return()
        }
    })
    
    # Run BRIC function
    bricOutput <- eventReactive(input$bric_launch, {
        if (input$bric_examplechoice == "bric_yes") {
            path <- "data/Yan_sub.txt"
        } else {
            path <- input$bric_file1$datapath
        }
        
        if (input$bric_method == "SC") {
            bric_K <- input$bric_celltype_number
            
        } else {
            bric_K <- NULL
        }
        
        # Crashes if application is restarted and console is not refreshed
        bric_out <- BRIC::final(
            i = path,
            method = input$bric_method,
            N = FALSE,
            R = FALSE,
            F = FALSE,
            d = FALSE,
            f = input$bric_f,
            k = input$bric_k,
            c = 1,
            o = input$bric_o,
            K = bric_K
        )

        return(
            data.frame(
                "cell_id" = names(bric_out),
                "cluster_number" = bric_out,
                row.names = seq_along(bric_out)
            )
        )
    })

    output$bric_output_table <- DT::renderDataTable({
        withProgress(mess = "Running BRIC...", value = 1, {
            bricOutput()
        })
    })
    
    ## DGE - Show download button (ALL)
    output$bric_download <- renderUI({
        validate(
            need(
                expr = !is.null(input$bric_launch),
                message = ""
            )
        )
        if(input$bric_launch == 0) {
            return()
        } else {
            downloadButton(
                "bric_download_data", "Download Cluster Predictions"
            )
        }
    })

    ## DGE - Download ALL data
    output$bric_download_data <- downloadHandler(
        filename = function() {
            paste("bric_cluster_pred_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            write.csv(
                bricOutput(),
                file, 
                row.names = FALSE
            )
        }
    ) 
}
