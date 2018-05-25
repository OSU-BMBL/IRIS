#---------------------------------------------------------------------
# Title:         IRIS - Shiny Application
# Author:        Brandon Monier
# Created:       2018-01-26 11:29:39 CDT
# Last Modified: 2018-05-25 11:31:14 CDT
#---------------------------------------------------------------------

## Set working directory (FOR LOCAL TESTING ONLY)
# setwd("D:/Box Sync/misc-github-repos/iris/")


# Load packages ----
source("pack-load.R")


# Server sources ----
source("iris-functions.R")
source("iris-xlsx.R")
source("tabs.R")


# Run shiny ----
shinyApp(irisUI, irisServer)
