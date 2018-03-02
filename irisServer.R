#---------------------------------------------------------------------
# Title:         IRIS - Server Script
# Author:        Brandon Monier
# Created:       2018-01-26 at 11:32:02
# Last Modified: 2018-01-26 at 11:32:39
#---------------------------------------------------------------------

# Change file upload size to 30 MB and sanitize errors
options(
  shiny.maxRequestSize = 30 * 1024^2,
  shiny.sanitize.errors = TRUE
)



# Server function
irisServer <- function(input, output) {
  
  ## Source exp. method functions
  source("iris-functions.R")

  ## Example data
  f1a <- as.matrix(read.csv("./data/count-data-small.csv", header=TRUE, row.names=1))
  f1b <- read.csv("./data/col-data-small.csv", header=TRUE, row.names=1)
  
  f2a <- as.matrix(read.csv("./data/count-data-big.csv", header=TRUE, row.names=1))
  f2b <- read.csv("./data/col-data-big.csv", header=TRUE, row.names=1)

  ## Data load option - count data
  output$file1 <- renderUI({
    if (input$examplechoice != "no") {
      return()
    } else if (input$examplechoice == "no") {
      fileInput(
        inputId = "file1", 
        label = "Submit count data (CSV)",
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      )      
    }
  })

  ## Data load option - metadata
  output$file2 <- renderUI({
    if (input$examplechoice != "no") {
      return()
    } else if (input$examplechoice == "no") {
      fileInput(
        inputId = "file2", 
        label = "Submit metadata (CSV)",
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      )      
    }
  })

	## QC - reactive - load and add data to DESeqDataSet class
	ddsout <- eventReactive(input$goqc, {
    if (input$examplechoice == "yes1") {
      cts <- f1a
      coldata <- f1b
    } else if (input$examplechoice == "yes2") {
      cts <- f2a
      coldata <- f2b
    } else if (input$examplechoice == "no") {
      cts <- input$file1
      coldata <- input$file2
      cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
      coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)      
    }

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

  ### QC - verbatim console - count data head
  output$fileoutputcts <- renderPrint(width = 400, {
    cts <- ddsout()[[4]]
    trubble(cts)
  })

  ### QC - verbatim console - metadata
  output$fileoutputcoldata <- renderPrint({
    coldata <- ddsout()[[2]]
    coldata
  })

  ### QC - verbatim console - gene no. (pre-filter)
  output$fileoutputcountpre <- renderPrint({
    cts.pre <- ddsout()[[4]]
    nrow(cts.pre)
  })

  ### QC - verbatim console - gene no. (pre-filter)
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


  ### QC - visualize - boxplot
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

  ### QC - visualize - histogram
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

  ### QC - visualize - barplot
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

  ### QC - select input - choose factor - PCA
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

  ### QC - select input - choose factor - MDS
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



  ##### B R E A K #####



  ## DEG - Choose analytical methods
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
        )
      )
    } else {
      selectInput(
        inputId = "dgemethod",
        label = "Choose method",
        choices = c(
          "DESeq2" = "deseq",
          "edgeR" = "edger"
        )
      )
    }
  })

  ## DEG - exp. setup 1 - two group comparisons - factor choice
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

  ## DEG - exp. setup 1 - two group comparisons - choose comparisons
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

  ## DEG - exp. setup 2 - mult. group comparisons - factor choice A
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

  ## DEG - exp. setup 2 - mult. group comparisons - factor choice B
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

  ## DEG - exp. setup 2 - mult. group comparisons - group comb. choices
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

  ## DEG - exp. setup 3 - interaction - factor choice A
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

  ## DEG - exp. setup 3 - interaction - factor choice B
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

  ## DEG - exp. setup 3 - interaction - reference level for factor A
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

  ## DEG - exp. setup 3 - interaction - reference level for factor B
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

  ## DEG - exp. setup 4 - additive model - blocking factor
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

  ## DEG - exp. setup 4 - additive model - treatment factor
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

  ## DEG - exp. setup 4 - additive model - reference level for blocking factor
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

  ## DEG - exp. setup 4 - additive model - reference level for treatment factor
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

  ## DEG - exp. setup 5 - main effect (ME) - choose ME
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

  ## DEG - exp. setup 5 - ME - choose ME reference
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

  ## DEG - exp. setup 5 - ME - choose contrasts
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

  ## DEG - exp. setup 6 - ME + group fact - choose ME
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

  ## DEG - exp. setup 6 - ME + group fact - choose group fact
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

  ## DEG - exp. setup 6 - ME + group fact - choose ME reference
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

  ## DEG - exp. setup 6 - ME + group fact - choose group fact level
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

  ## DEG - exp. setup 6 - ME group fact - choose contrasts
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

  ## DEG - exp. setup 7 - user input
  output$dgeexp7a <- renderUI({
    validate(
      need(input$dgemethod != "", "")
    ) 
    if (input$goqc == 0) {
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
  
  ## DEG - exp. setup - formula - header
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

  ## DEG - exp. setup - formula - formula 1
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

  ## DEG - exp. setup - formula - formula 2
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

  ## DEG - exp. setup 3 - interaction - formula layout
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

  ## DEG - exp. setup 4 - added effects - formula layout
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

  ## DEG - exp. setup 5 - ME - formula layout
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

  ## DEG - exp. setup 6A - ME + group fact. - formula layout
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

  ## DEG - exp. setup 6B - ME + group fact. - formula layout
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

  ## DEG - exp. setup 6C - ME + group fact. - formula layout
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

  ## DEG - exp. setup 6D - ME + group fact. - formula layout
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


  ## DEG - edgeR normalization option
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

  ## DEG - Choose contrasts
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

  ## DEG - user input for model matrix
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

  ## DEG - analysis - reactive
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

  ## DEG - select input - Visualization type
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
  
  ## DEG - Filter Data
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

  # DEG - Shared data
  share <- reactive({
    SharedData$new(dgeout3())
  })

  # DEG - Conditional plot
  output$dgeplot <- renderPlotly({
    validate(
      need(input$dgemaincontrasts != "", "")
    )
    check <- dgeout3()
    check <- nrow(check)
    if (input$godge == 0 | check < 1) {
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


  # DEG - Conditional table
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

  # DEG - Show download button (FILTERED)
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

  # Download FILTERED data
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
      write.csv(
        dgeout3()[input[["mytable_rows_all"]], ], 
        file, row.names = FALSE
      )
    }
  )   

  # DEG - Show download button (ALL)
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

  ## DEG - Download ALL data
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

  ## DEG - Generate overview data
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

  # DGE - Show download button (DGE Overview Table)
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
      ## DGE - Show download button - MA plot (PDF)
      if(input$goqc == 0) {
        return()
      } else {
        downloadButton("dldgemapdfimg", "Download MA Plot (PDF)")
      }     
    } else if (input$plottype == "volplot") {
      ## DGE - Show download button - Volcano plot (PDF)
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
      ## DGE - Show download button - MA plot (PNG)
      if(input$goqc == 0) {
        return()
      } else {
        downloadButton("dldgemapngimg", "Download MA Plot (PNG)")
      }     
    } else if (input$plottype == "volplot") {
      ## DGE - Show download button - Volcano plot (PNG)
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
  
  
  
  
  
  
  
  ### ブ レ ー ク  B R E A K  ブ レ ー ク ###



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

  ## QC - Show download button - Heatmap (PDF)
  output$dlqcheatplot1pdf <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      downloadButton("dlqcheatplot1pdfimg", "Download Static R Plot (PDF)")
    }     
  })
  
  ## QC - Download plot - Heatmap (PDF)
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
  
  ## QC - Show download button - Heatmap (PNG)
  output$dlqcheatplot1png <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      downloadButton("dlqcheatplot1pngimg", "Download Static R Plot (PNG)")
    }     
  })
  
  ## QC - Download plot - Heatmap (PNG)
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
  
  ### HEAT - select input - choose factor - HEATMAP
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

  ## QC - Show download button - Heat counts (PDF)
  output$dlqcheatplot2pdf <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      downloadButton("dlqcheatplot2pdfimg", "Download Static R Plot (PDF)")
    }     
  })
  
  ## QC - Download plot - Heat counts (PDF)
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
  
  ## QC - Show download button - Heat counts (PNG)
  output$dlqcheatplot2png <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      downloadButton("dlqcheatplot2pngimg", "Download Static R Plot (PNG)")
    }     
  })
  
  ## QC - Download plot - Heat counts (PNG)
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

  
  

  ### B R E A K ###



  
  
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
      actionButton("gobic", "Launch Analysis", icon = icon("space-shuttle"))
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

  ### BIC - head (3) - Clust overview
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
        " samples in cluster ", clust, ". To see what IDs were found in this cluster, click on the 'Download Cluster IDs' button at the bottom of the page." 
      )
      p(paste(out))
    }
  })

  ### BIC - head (3) - Clust overview
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


  ### B R E A K ###



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
        tooltips <- matrix(tooltips, ncol = ncol(cor.mat), byrow = FALSE)
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
      downloadButton("dlqccorplot1pdfimg", "Download Static R Plot (PDF)")
    }     
  })
  
  ## COR - DOWNLOAD PLOT - Corplot1 (pdf)
  output$dlqccorplot1pdfimg <- downloadHandler(
    filename =  function() {
      paste("cor-matrix.pdf")
    },
    content = function(file) {
      pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
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
      downloadButton("dlqccorplot1pngimg", "Download Static R Plot (PNG)")
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
      downloadButton("dlqccorplot2pdfimg", "Download Static R Plot (PDF)")
    }     
  })
  
  ## COR - DOWNLOAD PLOT - corplot2 (pdf)
  output$dlqccorplot2pdfimg <- downloadHandler(
    filename =  function() {
      paste("cor-scatterplot.pdf")
    },
    content = function(file) {
      pdf(file, width = 8, height = 6.5, onefile = FALSE) # open the pdf device
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
      downloadButton("dlqccorplot2pngimg", "Download Static R Plot (PNG)")
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
      pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
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



  ### B R E A K ###
  
  
  
  ## SESSION - Add session info verbatim
  output$sessinfo <- renderPrint({
    sessionInfo()
  })
}
