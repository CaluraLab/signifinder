% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureFunction.R
\name{IPSOVSign}
\alias{IPSOVSign}
\title{IPSOV Signature}
\usage{
IPSOVSign(
  dataset,
  nametype = "SYMBOL",
  inputType = "microarray",
  whichAssay = "norm_expr",
  ...
)
}
\arguments{
\item{dataset}{Normalized expression values. A data frame or a matrix
where rows correspond to genes and columns correspond to samples.
Alternatively, an object of type \linkS4class{SummarizedExperiment},
\code{\link[SingleCellExperiment]{SingleCellExperiment}} or
\code{\link[SpatialExperiment]{SpatialExperiment}} where the normalized
expression values should be in an assay called 'norm_expr'.}

\item{nametype}{character string saying the type of gene name ID (row names
in dataset). Either one of "SYMBOL", "ENTREZID" or "ENSEMBL".}

\item{inputType}{character string saying the type of data you are using.
Either one of "microarray" or "rnaseq".}

\item{whichAssay}{integer scalar or string indicating which assay of
dataset to use.}

\item{...}{other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
function.}
}
\value{
If dataset is a \linkS4class{SummarizedExperiment} object, then
scores are added in the \code{\link[SummarizedExperiment]{colData}} section.
If dataset is a data frame or a matrix, then a
\linkS4class{SummarizedExperiment} object is created in which scores are
added in the \code{\link[SummarizedExperiment]{colData}} section.
}
\description{
This signature is computed accordingly to the reference paper,
to have more details explore the function
\code{\link[signifinder]{availableSignatures}}.
}
\examples{
data(ovse)
IPSOVSign(dataset = ovse)

}
