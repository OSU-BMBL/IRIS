---
 Title:         IRIS Documentation
 Author:        Brandon Monier
 Created:       2018-01-12 at 09:56:43
 Last Modified: 2018-01-30 at 15:20:04
---

# IRIS - main effects experimental design

## About
This application is designed to test the functionality of main effects testing 
in addition to main effects with grouping factors.

## Get local app
Copy and paste this into your R environment:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "14-iris-me")
```