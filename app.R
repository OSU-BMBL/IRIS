#---------------------------------------------------------------------
# Title:         IRIS - Shiny Application
# Author:        Brandon Monier
# Created:       2018-01-26 11:29:39 CDT
# Last Modified: 2018-05-22 15:25:49 CDT
#---------------------------------------------------------------------

## Set working directory (FOR LOCAL TESTING ONLY)
# setwd("D:/Box Sync/misc-github-repos/iris/")


# Load packages ----
source("pack-load.R")


# Server sources ----
source("irisUI.R")
source("irisServer.R")
source("iris-xlsx.R")


# Run shiny ----
shinyApp(irisUI, irisServer)
