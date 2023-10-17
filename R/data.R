#' Example expression data.
#'
#' This is an example dataset containing gene expression values (in normalized
#' counts, TPM, CPM, and FPKM) of 40 ovarian cancer (OVC) patients extracted
#' from the Cancer Genome Atlas (TCGA) database.
#' This dataset should be used only with example purpose.
#' RNA sequencing OVC data were retrieved using
#' \code{\link[curatedTCGAData]{curatedTCGAData}} package. Data were then
#' normalized with the \code{\link[EDASeq]{betweenLaneNormalization}} function.
#' To lighten the dataset, the \code{\link[signifinder]{consensusOVSign}}
#' function was computed, which return 4 different scores, one for each OVC
#' subtype (Chen et al, 2018, Clinical Cancer Research) and the 10 samples
#' with the highest scores were selected for each subgroup.
#' Further, only the genes used for the signatures computation were kept.
#' Finally, all the signatures available in signifinder for OVC plus all the
#' pan-cancer signatures were computed.
#' Further details in signifinder/inst/scripts/howToGenerateOvse.Rmd.
#'
#' @docType data
#'
#' @return An object of class \linkS4class{SummarizedExperiment}.
#'
#' @usage data(ovse)
#'
"ovse"
