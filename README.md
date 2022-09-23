signifinder R package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

Signifinder is an R package that collects and implements 46
expression-based signatures from cancer literature. Through the analysis
of expression data with the collected signatures, signifinder can
attribute to each sample a score per signature that summarizes many
different tumor aspects, such as predict the response to therapy or the
survival association, as well as quantify multiple microenvironmental
conditions, such as hypoxia or the activity of the immune response.

## Installation

To install signifinder:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("signifinder")
```

<img src=./vignettes/figures/signifinder_visualization.jpg />
