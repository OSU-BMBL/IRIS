#---------------------------------------------------------------------
# Title:         IRIS - Tabs
# Author:        Brandon Monier
# Created:       2018-01-26 11:34:53 CDT
# Last Modified: 2018-05-25 14:54:13 CDT
#---------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
    title = "Welcome",
    icon = icon("hand-spock-o"),
    fluid = TRUE,
    sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
            width = 8,
            div(
                style = "display:inline-block;",
                style = "width:1000px;",
                style = "text-align:center;",
                style = "margin:auto;",
                tags$br(),
                tags$img(
                    src = "iris-logo-01.png",
                    width = "800px"
                ),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$button(
                    id = "iris_flow_01",
                    class = "btn action-button",
                    img(src = "iris-flow-01.png", height = "100px")
                ),
                tags$img(
                    src = "arrow-01.png",
                    height = "100px"
                ),
                tags$button(
                    id = "iris_flow_02",
                    class = "btn action-button",
                    img(src = "iris-flow-02.png", height = "100px")
                ),
                tags$br(),
                tags$img(
                    src = "arrow-02.png",
                    height = "20px"
                ),
                tags$br(),
                tags$button(
                    id = "iris_flow_03",
                    class = "btn action-button",
                    img(src = "iris-flow-03.png", height = "298px")
                ),
                tags$img(
                    src = "arrow-03.png",
                    height = "100px"
                ),
                tags$button(
                    id = "iris_flow_04",
                    class = "btn action-button",
                    img(src = "iris-flow-04.png", height = "298px")
                )
            ),
            tags$head(
                tags$style(
                    HTML(
                        "#iris_flow_01{background-color:white;}",
                        "#iris_flow_02{background-color:white;}",
                        "#iris_flow_03{background-color:white;}",
                        "#iris_flow_04{background-color:white;}"
                    )
                )
            ),
            bsTooltip(
                id = "iris_flow_01",
                title = paste(
                    "To use this server, submit",
                    "a raw read count matrix and a file for",
                    "sample treatment data. If you are performing single-",
                    "cell RNA-seq analysis, you will also need to provide a",
                    "file for ID lengths."
                ),
                placement = "top",
                options = list(container = "body")
            ),
            bsTooltip(
                id = "iris_flow_02",
                title = paste(
                    "Once you have submitted data, this section will aid",
                    "in the generation of a metadata file for GEO submission."
                ),
                placement = "top",
                options = list(container = "body")
            ),
            bsTooltip(
                id = "iris_flow_03",
                title = paste(
                    "You can also proceed here for initial analysis of your",
                    "count data."
                ),
                placement = "top",
                options = list(container = "body")
            ),
            bsTooltip(
                id = "iris_flow_04",
                title = paste(
                    "This section will provide the tools for differential",
                    "gene expression using 3 Bioconductor packages (DESeq2,",
                    "edgeR, or limma-voom) and also data export."
                ),
                placement = "top",
                options = list(container = "body")
            )
            # includeMarkdown("./markdown/welcome.md")
        )
    )
)



# Submit and QC page ----
tab.submit <- tabPanel(
    title = "Submit and QC ",
    icon = icon("filter"),
    fluid = TRUE,
    value = "val1",
    sidebarLayout(
        sidebarPanel(
            em(
                paste0(
                    "Confused? Consider reading the tutorial under: "
                )
            ),
            code("More -> Tutorial"),
            h4("1. Submission Parameters"),
            radioButtons(
                inputId = "examplechoice",
                label = "How do you want to start?",
                choices = c(
                    "Start with an example data set (small)." = "yes1",
                    "Start with an example data set (big)." = "yes2",
                    "Start with an example data set (scRNA)." = "yes3",
                    "Load my own data." = "no",
                    "Load my own data (scRNA)." = "scrna",
                    "Load my own data (scRNA - 10X Genomics)." = "scrna10x"
                )
            ),
            bsTooltip(
                id = "examplechoice",
                title = paste(
                    "Note: if you choose the big example data set or scRNA",
                    "example data set, analysis and visualization time will",
                    "be considerably longer than the small example data set."
                ),
                placement = "right",
                options = list(container = "body")
            ),
            uiOutput("file1"),
            uiOutput("file2"),
            uiOutput("file3"),
            br(),
            h4("2. Data Processing"),
            uiOutput("filter_choice"),
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
            actionButton(
                inputId = "goqc",
                "Submit",
                icon = icon("space-shuttle")
            ),
            br(),
            br(),
            p(
                paste(
                    "After you click \"submit\", you may proceed to either",
                    "the \"Data-Driven Analysis\", \"DGE Analysis\", or",
                    "\"GEO\" tabs for additional analyses."
                )
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
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqcboxplotpdf")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqcboxplotpng")
                    ),
                    br(),
                    br(),
                    br(),
                    uiOutput("counthist"),
                    plotlyOutput("hist"),
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqchistpdf")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqchistpng")
                    ),
                    br(),
                    br(),
                    br(),
                    uiOutput("counttotal"),
                    plotlyOutput("barplot"),
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqcbarplotpdf")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("dlqcbarplotpng")
                    ),
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
    title = "Discovery-Driven Analyses ",
    icon = icon("search"),
    fluid = TRUE,
    value = "val3",
    mainPanel(
        tabsetPanel(
            tabPanel(
                title = "Correlation",
                uiOutput("headcor"),
                plotlyOutput("corplot1"),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqccorplot1pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqccorplot1png")
                ),
                br(),
                br(),
                br(),
                plotlyOutput("corplot2"),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqccorplot2pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqccorplot2png")
                ),
                br(),
                br(),
                br(),
                uiOutput("headcor2"),
                plotOutput("corplot3", width = 600, height = 550),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcorplot3pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcorplot3png")
                ),
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
                title = "tSNE",
                uiOutput("headTSNE"),
                uiOutput("tsnefact"),
                uiOutput("tsneDim"),
                uiOutput("tsnePerp"),
                uiOutput("tsnePerpCheck"),
                plotlyOutput("tsnePlot"),
                br(),
                div(style = "display:inline-block", uiOutput("dlqctsnepdf")),
                div(style = "display:inline-block", uiOutput("dlqctsnepng")),
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
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcheatplot1pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcheatplot1png")
                ),
                br(),
                br(),
                br(),
                uiOutput("heatfactor"),
                plotlyOutput("heatplot2"),
                br(),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcheatplot2pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("dlqcheatplot2png")
                ),
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
                div(
                    style = "display:inline-block",
                    uiOutput("downloadbicfilt")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadbicplotpdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadbicplotpng")
                ),
                br(),
                br(),
                br()
            ),
            tabPanel(
                title = "Clustering",
                uiOutput("headclust"),
                uiOutput("headclustwarn"),
                uiOutput("clustvarnumber"),
                uiOutput("clustalg"),
                uiOutput("goclust"),
                br(),
                br(),
                uiOutput("headclustplotW01"),
                plotOutput("clustplotW01"),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW01png")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW01pdf")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustmodK")
                ),
                uiOutput("headclustplotW02"),
                plotOutput("clustplotW02"),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW02png")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW02pdf")
                ),
                uiOutput("headclustplotW03"),
                plotOutput("clustplotW03"),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW03png")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustplotW03pdf")
                ),
                br(),
                br(),
                uiOutput("headclustmoddown"),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustmod")
                ),
                div(
                    style = "display:inline-block",
                    uiOutput("downloadclustsample")
                ),
                br(),
                br(),
                br()
            )
        )
    )
)


# DEG Analysis ----
tab.deg <- tabPanel(
    title = "DGE Analysis ",
    icon = icon("bar-chart"),
    fluid = TRUE,
    value = "val4",
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
                    div(
                        style = "display:inline-block",
                        uiOutput("dldgeoverpdf")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("dldgeoverpng")
                    ),
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
                    div(
                        style = "display:inline-block",
                        uiOutput("dldgemavolpdf")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("dldgemavolpng")
                    ),
                    br(),
                    br(),
                    br(),
                    DT::dataTableOutput("mytable"),
                    br(),
                    div(
                        style = "display:inline-block",
                        uiOutput("downloadfilt")
                    ),
                    div(
                        style = "display:inline-block",
                        uiOutput("downloadall")
                    ),
                    br(),
                    br(),
                    br()
                ),
                tabPanel(
                    title = "Functional Enrichment",
                    includeMarkdown("./markdown/functional-enrichment.md")
                )
            )
        )
    )
)



# GEO ----
tab.geo <- tabPanel(
    title = "GEO",
    icon = icon("database"),
    fluid = TRUE,
    value = "val2",
    sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
            includeMarkdown("./markdown/geo.md"),
            br(),
            h3("Questionnaire"),
            uiOutput("geo_series"),
            uiOutput("geo_01_title"),
            uiOutput("geo_02_summary"),
            uiOutput("geo_xx_overall_design"),
            div(id = "contrib"),
            uiOutput("geo_03a_contrib_action"),
            uiOutput("geo_03b_contrib"),
            uiOutput("geo_04_suppfile"),
            uiOutput("geo_05_sra"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_series")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_series_check")
            ),
            uiOutput("geo_space_01"),


            uiOutput("geo_samples"),
            uiOutput("geo_samples_list"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_sample")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_sample_check")
            ),
            uiOutput("geo_space_02"),


            uiOutput("geo_protocols"),
            uiOutput("geo_07_growth_protocol"),
            uiOutput("geo_08_treatment_protocol"),
            uiOutput("geo_09_extract_protocol"),
            uiOutput("geo_10_lib_construct_protocol"),
            uiOutput("geo_11_lib_strategy"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_protocol")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_protocol_check")
            ),
            uiOutput("geo_space_03"),


            uiOutput("geo_data_proc_pipeline"),
            div(id = "data_step"),
            uiOutput("geo_12a_data_proc_action"),
            uiOutput("geo_13_genome_build"),
            uiOutput("geo_14_proc_data_files"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_pipeline")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_pipeline_check")
            ),
            uiOutput("geo_space_04"),


            uiOutput("geo_proc_head"),
            uiOutput("geo_proc_data"),
            uiOutput("geo_proc_data_list"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_proc_data")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_proc_data_check")
            ),
            uiOutput("geo_space_05"),


            uiOutput("geo_raw_head"),
            uiOutput("geo_raw_data"),
            uiOutput("geo_raw_data_list"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_raw_data")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_raw_data_check")
            ),
            uiOutput("geo_space_06"),


            uiOutput("geo_pair_head"),
            uiOutput("geo_paired_end"),
            uiOutput("geo_paired_end_list"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_paired_end")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_paired_end_check")
            ),
            uiOutput("geo_space_07"),


            uiOutput("geo_solid_head"),
            uiOutput("geo_solid"),
            uiOutput("geo_solid_list"),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_solid")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_submit_solid_check")
            ),
            uiOutput("geo_space_08"),



            uiOutput("geo_downloads"),
            div(
                style = "display:inline-block",
                uiOutput("geo_download_excel_1")
            ),
            div(
                style = "display:inline-block",
                uiOutput("geo_download_proc_data_1")
            ),
            # verbatimTextOutput("geo_debug"),
            br(),
            br(),
            br(),
            br(),
            br()
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
