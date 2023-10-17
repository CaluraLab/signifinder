\name{NEWS}
\title{ for package \pkg{signifinder} }

\section{Version 1.4.0, Bioconductor 3.18 Release}{
 \itemize{
  \item Add \code{glioCellStateSign} function to compute the glioblastoma
  cellular states defined by Neftel C. et al. Cell (2019).
  \item Add \code{whichAssay} argument to signature functions to allow the
  user to specify which assay to use for the signature computation.
 }
}

\section{Version 1.2.1, Bioconductor 3.17 Release}{
 \itemize{
  \item Add \code{evaluationSignPlot} function to show some technical
  information of the signatures computed.
  \item Add \code{nametype} argument to \code{geneHeatmapSignPlot} function
  to allow more gene name ID in data.
  \item The vignette now contains an example with a single-cell dataset and
  an example with a spatial transcriptomics dataset.
 }
}

\section{Version 1.0.0, Bioconductor 3.16 Release}{
 \itemize{
  \item New package signifinder, for implementation of transcriptional cancer
  signatures.
 }
}
