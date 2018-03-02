#---------------------------------------------------------------------
# Title:         IRIS - Tabs
# Author:        Brandon Monier
# Created:       2018-01-26 at 11:34:53
# Last Modified: 2018-01-26 at 11:35:57
#---------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
  title = "Welcome", icon = icon("hand-spock-o"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      includeMarkdown("./markdown/welcome.md")
    )
  )
)



# Submit and QC page ----
tab.submit <- tabPanel(
  title = "Submit and QC ", icon = icon("filter"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("1. Submission Parameters"),
      radioButtons(
        inputId = "examplechoice",
        label = "How do you want to start?",
        choices = c(
          "Start with a small example data set." = "yes1",
          "Start with a big example data set." = "yes2",
          "Load my own data." = "no"  
        )
      ),
      bsTooltip(
        id = "examplechoice", 
        title = "Note: if you choose the big example data set, analysis and visualization time will be considerably longer than the small example data set.",
        placement = "bottom"
      ),
      uiOutput("file1"),
      uiOutput("file2"),
      br(),
      h4("2. Data Processing"),
      textInput(
        inputId = "prefilt",
        label = withMathJax("Filter cutoff (count data row sums \\( < n\\))"),
        value = as.numeric(10)
      ),
      selectInput(
        inputId = "transform",
        label = "Choose transformation method for counts",
        choices = c(
          "Normal log: log2(n + pseudocount)" = "log",
          "Regularized log: rlog(n)" = "rlog",
          "Variance stabilizing transform: vst(n)" = "vst",
          "No transformation" = "raw"
        )       
      ),
      br(),
      h4("3. Launch Overview"),
      actionButton(inputId = "goqc", "Submit", icon = icon("space-shuttle")),
      br(),
      br(),
      p(
        "After you click 'submit', you may proceed to either the 'Preliminary Analysis' or 'DGE Analysis' tab for additional analyses."
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "File Summary",
          uiOutput("filesummarycts"),
          verbatimTextOutput("fileoutputcts"),
          uiOutput("filesummarycoldata"),
          verbatimTextOutput("fileoutputcoldata"),
          br(),
          uiOutput("headcountpre"),
          verbatimTextOutput("fileoutputcountpre"),
          uiOutput("headcountpost"),
          verbatimTextOutput("fileoutputcountpost"),
          br(),
          br()
        ),
        tabPanel(
          title = "Count Summary",
          uiOutput("countbox"),
          plotlyOutput("boxplot"),
          div(style = "display:inline-block", uiOutput("dlqcboxplotpdf")),
          div(style = "display:inline-block", uiOutput("dlqcboxplotpng")),
          br(),
          br(),
          br(),
          uiOutput("counthist"),
          plotlyOutput("hist"),
          div(style = "display:inline-block", uiOutput("dlqchistpdf")),
          div(style = "display:inline-block", uiOutput("dlqchistpng")),
          br(),
          br(),
          br(),
          uiOutput("counttotal"),
          plotlyOutput("barplot"),
          div(style = "display:inline-block", uiOutput("dlqcbarplotpdf")),
          div(style = "display:inline-block", uiOutput("dlqcbarplotpng")),
          br(),
          br(),
          br()
        )
      )
    )
  )
)



# Exploratory Analysis
tab.prelim <- tabPanel(
  title = "Preliminary Analysis ", icon = icon("search"),
  fluid = TRUE,
  mainPanel(
    tabsetPanel(
      tabPanel(
        title = "Correlation",
        uiOutput("headcor"),
        plotlyOutput("corplot1"),
        div(style = "display:inline-block", uiOutput("dlqccorplot1pdf")),
        div(style = "display:inline-block", uiOutput("dlqccorplot1png")),
        br(),
        br(),
        br(),
        plotlyOutput("corplot2"),
        div(style = "display:inline-block", uiOutput("dlqccorplot2pdf")),
        div(style = "display:inline-block", uiOutput("dlqccorplot2png")),
        br(),
        br(),
        br(),
        uiOutput("headcor2"),
        plotOutput("corplot3", width = 600, height = 550),
        div(style = "display:inline-block", uiOutput("dlqcorplot3pdf")),
        div(style = "display:inline-block", uiOutput("dlqcorplot3png")),
        br(),
        br(),
        br()
      ),
      tabPanel(
        title = "PCA",
        uiOutput("headpca"),
        uiOutput("pcafact"),
        plotlyOutput("pca"),
        br(),
        div(style = "display:inline-block", uiOutput("dlqcpcapdf")),
        div(style = "display:inline-block", uiOutput("dlqcpcapng")),
        br(),
        br(),
        br()
      ),
      tabPanel(
        title = "MDS",
        uiOutput("headmds"),
        uiOutput("mdsfact"),
        plotlyOutput("mds"),
        br(),
        div(style = "display:inline-block", uiOutput("dlqcmdspdf")),
        div(style = "display:inline-block", uiOutput("dlqcmdspng")),
        br(),
        br(),
        br()
      ),
      tabPanel(
        title = "Heatmap",
        uiOutput("headheat"),
        uiOutput("heatnumber"),
        plotlyOutput("heatplot1"),
        br(),
        div(style = "display:inline-block", uiOutput("dlqcheatplot1pdf")),
        div(style = "display:inline-block", uiOutput("dlqcheatplot1png")),
        br(),
        br(),
        br(),
        uiOutput("heatfactor"),
        plotlyOutput("heatplot2"),
        br(),
        div(style = "display:inline-block", uiOutput("dlqcheatplot2pdf")),
        div(style = "display:inline-block", uiOutput("dlqcheatplot2png")),
        br(),
        br(),
        br()
      ),
      tabPanel(
        title = "Biclustering",
        uiOutput("headbic"),
        uiOutput("headbicparameters"),
        uiOutput("bicvarnumber"),
        uiOutput("bicalg"),
        uiOutput("gobic"),
        br(),
        br(),
        uiOutput("headbicsummary"),
        uiOutput("bicclustnumber"),
        uiOutput("bicsummary"),
        br(),
        uiOutput("headbicheatsummary"),
        plotOutput("bicheatplot", height = 800, width = 600),
        div(style = "display:inline-block", uiOutput("downloadbicfilt")),
        div(style = "display:inline-block", uiOutput("downloadbicplotpdf")),
        div(style = "display:inline-block", uiOutput("downloadbicplotpng")),
        br(),
        br(),
        br()
      )
    )
  )
)


# DEG Analysis ----
tab.deg <- tabPanel(
  title = "DGE Analysis ", icon = icon("bar-chart"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("1. Experimental Setup"),
      selectInput(
        inputId = "dgeexpsetup",
        label = "Choose an experimental design",
        choices = c(
          "Two group comparisons" = "exp1",
          "Multiple factor comparisons (factorial)" = "exp2",
          "Classical interaction design" = "exp3",
          "Additive models - paired or blocking" = "exp4",
          "Main effects" = "exp5",
          "Main effects with grouping factors" = "exp6",
          "Custom design (advanced)" = "exp7"
        ),
        selected = "",
        multiple = FALSE
      ),
      br(),
      h4("2. DGE Parameters"),
      uiOutput("dgemethod"),
      splitLayout(
        textInput(
          inputId = "dgepadjcutoff",
          label = withMathJax("Adj. \\(p\\)-value cutoff"),
          value = 0.05
        ),
        textInput(
          inputId = "dgefcmin",
          label = "Min. fold change",
          value = 1
        )
      ),
      uiOutput("dgeexpedgernorm"),
      br(),
      h4("3. Experimental Parameters"),
      uiOutput("dgeexp1a"),
      uiOutput("dgeexp1b"),
      uiOutput("dgeexp2a"),
      uiOutput("dgeexp2b"),
      uiOutput("dgeexp2c"),
      uiOutput("dgeexp3a"),
      uiOutput("dgeexp3b"),
      uiOutput("dgeexp3c"),
      uiOutput("dgeexp3d"),
      uiOutput("dgeexp4a"),
      uiOutput("dgeexp4b"),
      uiOutput("dgeexp4c"),
      uiOutput("dgeexp4d"),
      uiOutput("dgeexp5a"),
      uiOutput("dgeexp5b"),
      uiOutput("dgeexp5c"),
      uiOutput("dgeexp6a"),
      uiOutput("dgeexp6b"),
      uiOutput("dgeexp6c"),
      uiOutput("dgeexp6d"),
      uiOutput("dgeexp6e"),
      uiOutput("dgeexp7a"),
      uiOutput("dgeexpformhead"),
      uiOutput("dgeexpform1"),
      uiOutput("dgeexpform2"),
      uiOutput("dgeexpform3"),
      uiOutput("dgeexpform4"),
      uiOutput("dgeexpform5"),
      uiOutput("dgeexpform6a"),
      uiOutput("dgeexpform6b"),
      uiOutput("dgeexpform6c"),
      uiOutput("dgeexpform6d"),
      br(),
      h4("4. Launch Analysis"),
      actionButton("godge", "Submit", icon = icon("space-shuttle"))
    ),
    mainPanel = mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Overview",
          uiOutput("headdgeoverview"),
          verbatimTextOutput("debugdge2"),
          uiOutput("dgeoverview"),
          uiOutput("dldgeoverviewtbl"),
          br(),
          br(),
          br(),
          plotlyOutput("dgeplot2", height = 600),
          br(),
          div(style = "display:inline-block", uiOutput("dldgeoverpdf")),
          div(style = "display:inline-block", uiOutput("dldgeoverpng")),
          br(),
          br(),
          br()
        ),
        tabPanel(
          title = "Plots",
          uiOutput("headdgeplots"),
          uiOutput("dgemaincontrasts"),
          verbatimTextOutput("debugdge"),
          uiOutput("vistype"),
          plotlyOutput("dgeplot"),
          br(),
          div(style = "display:inline-block", uiOutput("dldgemavolpdf")),
          div(style = "display:inline-block", uiOutput("dldgemavolpng")),
          br(),
          br(),
          br(),
          DT::dataTableOutput("mytable"),
          br(),
          div(style = "display:inline-block", uiOutput("downloadfilt")),
          div(style = "display:inline-block", uiOutput("downloadall")),
          br(),
          br(),
          br()
        )
      )
    )
  )
)


# Help ----
tab.tutorial <- tabPanel(
  title = "Tutorial",
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      includeMarkdown("./markdown/iris-tutorial.md")
    )
  )
)


# FAQ ----
tab.faq <- tabPanel(
  title = "FAQ",
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      includeMarkdown("./markdown/faq.md")
    )
  )
)


# About Us ----
tab.about <- tabPanel(
  title = "About Us",
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      includeMarkdown("./markdown/about-us.md")
    )
  )
)


# System Info ----
tab.sessinfo <- tabPanel(
  title = "Session Info",
  flud = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      h2("R Session Info"),
      verbatimTextOutput("sessinfo")
    )
  )
)


