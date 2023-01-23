signifinder R package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

`signifinder` is an R package that collects and implements 46
expression-based signatures from cancer literature. Through the analysis
of expression data with the collected signatures, `signifinder` can
attribute to each sample a score per signature that summarizes many
different tumor aspects, such as predict the response to therapy or the
survival association, as well as quantify multiple microenvironmental
conditions, such as hypoxia or the activity of the immune response.

The `signifinder` package is available from
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/signifinder.html),
where a vignette containing examples and documentation is included.

## Installation

The `signifinder` package can be installed from Bioconductor:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("signifinder")
```

<img src=./vignettes/figures/signifinder_main_figure.png />
