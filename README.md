
<!-- README.md is generated from README.Rmd. Please edit that file -->
hbvPSO (R package)
==================

hbvPSO provides an easy way to apply particle swarm optimization to the HBV rainfall-runoff model. It uses the [`hydroPSO`](https://cran.r-project.org/web/packages/hydroPSO/index.html) package for the optimization and the package [`TUWmodel`](https://cran.r-project.org/web/packages/TUWmodel/index.html) for the model itself. Input parameters can optionally be supplied as configuration files written in R, and multiple optimization runs can be executed in batches. Contains some helper functions to import data and settings from HBV-light.

Installation
------------

You can install hbvPSO from github with:

``` r
# install.packages("devtools")
devtools::install_github("jthurner/hbvPSO")
```

Status
------

This package is still work-in-progress. Using multiple zones does not work reliably at the moment. You might also need to install patched versions of some dependencies to get full functionality:

``` r
# required to allow plotting of sim-obs graph:
devtools::install_github("jthurner/hydroPSO",ref = "plotting")
# adds lnNSE to the objective functions in the sim-obs graph: 
devtools::install_github("hzambran/hydroGOF",ref = "ggof-lnnse")
```

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- ## basic example code -->
<!-- ``` -->
