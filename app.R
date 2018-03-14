#---------------------------------------------------------------------
# Title:         IRIS - Shiny Application
# Author:        Brandon Monier
# Created:       2018-01-26 at 11:29:39
# Last Modified: 2018-01-26 at 11:31:06
#---------------------------------------------------------------------

## Set working directory (FOR LOCAL TESTING ONLY)
setwd("D:/Box Sync/misc-github-repos/iris/")


# Load packages ----
source("pack-load.R")


# Server sources ----
source("irisUI.R")
source("irisServer.R")


# Run shiny ----
shinyApp(irisUI, irisServer)
