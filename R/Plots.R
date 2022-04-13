
#' Scatterplot for a single signature
#'
#' Given signatures' scores it returns a scatterplot of patients scores and a
#' barplot of the density distributions of patients scores.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions.
#' @param whichSign The signature function to plot.
#' @param statistics The statistics to be plot in the graph. It should be a
#' pre-defined character between "mean", "median" or "quantiles".
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median quantile
#'
#' @export
oneSignPlot <- function(data, whichSign, statistics = NULL){

    if(length(whichSign)>1){
        stop("you must provide only one signature for this plot")}

    signatureNameCheck(data, whichSign)

    if(!(is.null(statistics))){
        if(!(statistics %in% c("mean", "median", "quantiles"))){
            stop("statistics must be one of: mean, median and quantiles.")}}

    signval <- sort(data.frame(colData(data))[,whichSign])

    g1 <- ggplot() +
        geom_point(mapping = aes(signval, seq_along(signval)), size = 1) +
        labs(x = whichSign, y = "Index") +
        theme(panel.background = element_blank(),
              axis.line = element_line(colour = "grey50"))
    g2 <- ggplot() +
        geom_histogram(mapping = aes(signval, ..density..),
                       colour = "white", fill = "skyblue2") +
        geom_density(mapping = aes(signval), size = 1.5) +
        labs(x = whichSign, y = "Density") +
        theme(panel.background = element_blank(),
              axis.line = element_line(colour = "grey50"))

    if(!is.null(statistics)){
        if(statistics=="mean"){
            g1 <- g1 + geom_vline(mapping = aes(xintercept = mean(signval)),
                                  col = "red", size = 1.2, linetype = 2)
            g2 <- g2 + geom_vline(mapping = aes(xintercept = mean(signval)),
                                  col = "red", size = 1.2, linetype = 2)}
        else if (statistics=="median"){
            g1 <- g1 + geom_vline(mapping = aes(xintercept = median(signval)),
                                  col = "red", size = 1.2, linetype = 2)
            g2 <- g2 + geom_vline(mapping = aes(xintercept = median(signval)),
                                  col = "red", size = 1.2, linetype = 2)}
        else if (statistics=="quantiles"){
            g1 <- g1 +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.25)),
                           col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.50)),
                           col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.75)),
                           col = "red", size = 1.2, linetype = 2)
            g2 <- g2 +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.25)),
                           col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.50)),
                           col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.75)),
                           col = "red", size = 1.2, linetype = 2)}}
    g1 <- g1 +
        annotate("text", label = statistics, x = quantile(signval, 0.10),
                 y = length(signval), size = 4, colour = "red")

    return(g1 + g2)
}


#' Gene's Signatures' Heatmap
#'
#' Given one or multiple signatures, the function gives the opportunity to
#' observe the trend of the signature's scores based on the expression of the
#' genes included in each of them.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions.
#' @param whichSign The signatures to show in the plot. This will automatically
#' select corresponding genes to show in the plot.
#' @param logCount Logical. Give the opportunity to compute logarithm for the
#' expression values in the dataset.
#' @param splitBySign Logical. Give you the possibility to split heatmap's row
#' based on the signature's origin of the genes.
#' @param sampleAnnot A vector containing samples annotations.
#' @param splitByAnnot A categorical variable. Give you the possibility to split
#' columns based on samples' annotations.
#' @param ... Other parameters specific of the function
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation '%v%' HeatmapAnnotation
#' @importFrom magrittr '%>%'
#' @importFrom dplyr group_by summarise_all
#' @importFrom SummarizedExperiment colData
#'
#' @export
geneHeatmapSignPlot <- function(data, whichSign, logCount = FALSE,
                                splitBySign = FALSE, sampleAnnot = NULL,
                                splitByAnnot = FALSE, ...){

    signatureNameCheck(data, whichSign)

    dataset <- getMatrix(data)

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=ncol(dataset)){
            stop("sampleAnnot length is different than samples dimension")}
    } else {if(splitByAnnot){
        stop("splitByAnnot can be TRUE only if sampleAnnot is provided")}}

    signval <- colData(data)[,whichSign]
    if(length(whichSign)==1){
        signval <- matrix(
            signval, nrow = 1, dimnames = list(whichSign, colnames(dataset)))
    } else {
        signval <- sapply(signval, range01)
        signval <- as.matrix(t(signval))}

    geneTable <- as.data.frame(do.call(rbind, lapply(whichSign, GetGenes))) %>%
        group_by(Gene) %>% summarise_all(paste, collapse=",")
    signatureGenes <- geneTable$Gene

    filtdataset <- as.matrix(dataset[row.names(dataset) %in% signatureGenes, ])

    dots <- list(...)
    htargs <- matchArguments(dots, list(
        name = "Genes", show_column_names = FALSE, col = mycol,
        row_names_gp = grid::gpar(fontsize = 6)))

    if(logCount){htargs$matrix = log2(filtdataset+1)
    } else {htargs$matrix = filtdataset}

    if(length(whichSign)!=1){
        signAnnot <- geneTable$Signature[geneTable$Gene %in% rownames(filtdataset)]
        if(splitBySign){
            htargs$row_split = signAnnot
        } else {
            ha <- rowAnnotation(signature = signAnnot)
            htargs$right_annotation = ha}}

    if(splitByAnnot & is.character(sampleAnnot)){
        htargs$column_split = sampleAnnot
        ht <- Heatmap(signval, name = "Signature", col = mycol1,
                      column_split = sampleAnnot)
    } else {
        if(!is.null(sampleAnnot)){
            hatop = HeatmapAnnotation(sampleAnnot = sampleAnnot)
            htargs$top_annotation = hatop}
        ht <- Heatmap(signval, name = "Signature", col = mycol1)}

    ht2 <- do.call(Heatmap, htargs)
    g <- ht %v% ht2

    return(g)
}


#' Global Heatmap of Signatures' scores.
#'
#' Given a matrix of multiple signatures, the function gives the opportunity to
#' observe the trend of the signature's scores for each samples. Moreover it
#' gives the possibility to order all the matrix based on a single signature.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions. Th function takes as input a matrix where each
#' row represents a signature and each column represent a samples.
#' @param whichSign The signatures to show in the plot.
#' @param clusterBySign One ore more signatures that clusterize columns on
#' Heatmap.
#' @param sampleAnnot A vector containing samples annotations.
#' @param splitByAnnot A categorical variable. Give you the possibility to split
#' columns based on samples' annotations.
#' @param ... Other parameters specific of the function
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @importFrom ComplexHeatmap Heatmap '%v%'
#' @importFrom SummarizedExperiment colData
#'
#' @export
heatmapSignPlot <- function(data, whichSign = NULL, clusterBySign = NULL,
                            sampleAnnot = NULL, splitByAnnot = FALSE, ...){

    if(!is.null(whichSign)){signatureNameCheck(data, whichSign)}
    if(!is.null(clusterBySign)){signatureNameCheck(data, clusterBySign)}

    if(sum(colnames(colData(data)) %in% SignatureNames)>0){
        data <- colData(data)[, colnames(colData(data)) %in% SignatureNames]
    } else {stop("There are no signatures in data")}

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=nrow(data)){
            stop("sampleAnnot length is different than samples dimension")}
    } else {if(splitByAnnot){
        stop("splitByAnnot can be TRUE only if sampleAnnot is provided")}}

    if(!is.null(whichSign)){
        data <- data[, intersect(colnames(data), c(whichSign, clusterBySign))]}
    keepnames <- rownames(data)

    data <- sapply(data, range01)
    row.names(data) <- keepnames
    data <- as.matrix(t(data))

    dots <- list(...)
    htargs <- matchArguments(
        dots, list(name = "Signatures", show_column_names = FALSE, col = mycol))

    if(!is.null(sampleAnnot)){
        if(splitByAnnot){
            htargs$column_split = sampleAnnot
        } else {
            hatop = HeatmapAnnotation(sampleAnnot = sampleAnnot)
            htargs$top_annotation = hatop}}

    if(is.null(clusterBySign)){
        htargs$matrix = data
        g <- do.call(Heatmap, htargs)
    } else {
        n <- which(rownames(data) %in% clusterBySign)
        fm <- as.matrix(data.frame(data)[n,])
        sm <- as.matrix(data.frame(data)[-n,])
        htargs$matrix = sm
        if(splitByAnnot){
            ht <- Heatmap(fm, name = "Guiding Signatures", col = mycol1,
                          column_split = sampleAnnot)
        } else {
            ht <- Heatmap(fm, name = "Guiding Signatures", col = mycol1)}
        g <- ht %v% do.call(Heatmap, htargs)
    }
    return(g)
}


#' Correlation Plot
#'
#' Given a matrix of multiple signatures, the function gives the opportunity to
#' observe signatures correlates between each other.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions.
#' @param whichSign The signatures to show in the plot.
#' @param sampleAnnot A categorical variable containing samples annotations.
#' @param selectByAnnot A subgroup from `sampleAnnot` to use to construct the
#' correlation plot
#'
#' @return A correlation ellipse graph.
#'
#' @importFrom openair corPlot
#' @importFrom SummarizedExperiment colData
#'
#' @export
correlationSignPlot <- function(data, whichSign = NULL, sampleAnnot = NULL,
                                selectByAnnot = NULL){

    if(!is.null(whichSign)){signatureNameCheck(data, whichSign)}

    tmp <- colData(data)

    if(sum(colnames(tmp) %in% SignatureNames)>0){
        if(is.null(whichSign)){
            signs <- intersect(SignatureNames, colnames(tmp))
        } else {
            signs <- Reduce(
                intersect, list(SignatureNames, colnames(tmp), whichSign))}
    } else {stop("There are no signatures in data")}

    tmp <- tmp[,signs]

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=nrow(tmp)){
            stop("sampleAnnot length is different than samples dimension")}
        if(!is.null(selectByAnnot)){
            if(!(selectByAnnot %in% sampleAnnot)){
                stop("selectByAnnot is not present in sampleAnnot")}
        } else {
            stop("sampleAnnot can be used only if selectByAnnot is also provided")}
    } else {if(!is.null(selectByAnnot)){
        stop("selectByAnnot can be used only if sampleAnnot is also provided")}}

    if(!is.null(sampleAnnot)){
        if(!is.null(selectByAnnot)){tmp <- tmp[sampleAnnot==selectByAnnot,]}}

    SignMatrix <- sapply(tmp, range01)

    g <- corPlot(as.data.frame(SignMatrix), cluster = TRUE,
                 dendrogram = TRUE, lower = TRUE, fontsize=7)

    return(g)
}


#' Survival Plot
#'
#' It creates survival curves from either a formula (e.g. the
#' Kaplan-Meier), a previously fitted Cox model, or a previously fitted
#' accelerated failure time model.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions.
#' @param survData A dataframe where each row represents a samples (equal to
#' the names in the expression matrix) and two columns: the first holding
#' survival data of time, indicating the follow up time of the samples; the
#' second holding data of the survival status, an indicator normally 0=alive and
#' 1=dead. For interval censored data, the status indicator is 0=right censored,
#' 1=event at time, 2=left censored, 3=interval censored.
#' @param whichSign The signatures to test in the plot.
#' @param cutpoint A character between "median", "mean" and "optimal" or a
#' numeric value, which divide samples between high values and low values. The
#' function computes the value with the method indicated or employs the values
#' directly supplied by the user, and based on that number it divides samples in
#' higher and lower to compare them as the two groups of the survival plot.
#' @param sampleAnnot A categorical variable containing samples annotations
#' named with samples names equal to the row names used in `survData`.
#' @param selectByAnnot A group from `sampleAnnot` used to compute the
#' survival analysis.
#'
#' @return A Survival plot and the statistics computed on data.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median
#'
#' @export
survivalSignPlot <- function(data, survData, whichSign, cutpoint = "mean",
                             sampleAnnot = NULL, selectByAnnot = NULL){

    if(length(whichSign)>1){
        stop("you must provide only one signature for this plot")}
    signatureNameCheck(data, whichSign)

    if(ncol(survData)>2){
        stop("survData should contain only two columns: survival data and status")}

    if(!(is.numeric(cutpoint))){
        if(!(cutpoint %in% c("mean", "median", "optimal"))){
            stop("Cutpoint must be either a number or one of: mean, median and optimal")}}

    tmp <- intersect(rownames(colData(data)), rownames(survData))
    tmp <- as.data.frame(cbind(colData(data)[tmp,], survData[tmp,]))

    grp  <- rep("high", nrow(tmp))
    names(grp) <- rownames(tmp)
    n <- ncol(tmp)
    colnames(tmp)[c(n-1,n)] <- c("survival", "status")
    if(cutpoint=="mean"){
        grp[which(tmp[,whichSign] < mean(tmp[,whichSign]))] <- "low"
    } else if(cutpoint=="median"){
        grp[which(tmp[,whichSign] < median(tmp[,whichSign]))] <- "low"
    } else if(cutpoint=="optimal"){
        optval <- maxstat::maxstat.test(
            survival::Surv(survival, status) ~ tmp[,whichSign],
            data = tmp, smethod="LogRank", pmethod="Lau94")
        grp[which(tmp[,whichSign] < optval$estimate)] <- "low"
    } else {grp[which(tmp[,whichSign] < cutpoint)] <- "low"}

    if((sum(grp=="low")<length(grp)/10) | (sum(grp=="low")>length(grp)*9/10)){
        warning(paste("groups size is non homogeneous:",
                      sum(grp=="low"), "low and",
                      sum(grp=="high"), "high"))}

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=nrow(tmp)){
            stop("sampleAnnot length is different than samples dimension")}
        if(!is.null(selectByAnnot)){
            if(!(selectByAnnot %in% sampleAnnot)){
                stop("selectByAnnot is not present in sampleAnnot")}
        } else {
            stop("sampleAnnot can be used only if selectByAnnot is also provided")}
    } else {if(!is.null(selectByAnnot)){
        stop("selectByAnnot can be used only if sampleAnnot is also provided")}}

    tmp <- cbind(tmp, grp)
    if(!is.null(sampleAnnot)){
        if(!is.null(selectByAnnot)){tmp <- tmp[sampleAnnot==selectByAnnot,]}}

    fit <- survival::survfit(survival::Surv(survival, status) ~ grp, data = tmp)

    g <- survminer::ggsurvplot(
        fit, data = tmp, risk.table = TRUE, legend.title = whichSign,
        palette = c("red", "blue"), ggtheme = ggplot2::theme_gray(15),
        font.legend = 15, font.tickslab = 15, font.x = 15, font.y = 15,
        risk.table.fontsize = 5, pval = TRUE, surv.median.line = "hv",
        risk.table.col = "strata", tables.height = 0.4)
    return(g)
}


#' Ridgeline Plot
#'
#' The function calculates densities from the point data mapped onto the x axis,
#' arranging multiple density plots of the Signatures' scores.
#'
#' @param data An object of type \linkS4class{SummarizedExperiment} coming from
#' signatures functions.
#' @param whichSign The signatures to test in the ridgeline plot.
#' @param groupByAnnot  A categorical variable containing samples annotations.
#' @param selectByAnnot A character indicating the group/s defined in
#' `groupByAnnot` to show in the plot.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom SummarizedExperiment colData
#'
#' @export
ridgelineSignPlot <- function(data, whichSign = NULL, groupByAnnot = NULL,
                              selectByAnnot = NULL){

    if(!is.null(whichSign)){signatureNameCheck(data, whichSign)}

    tmp <- colData(data)

    if(sum(colnames(tmp) %in% SignatureNames)>0){
        if(is.null(whichSign)){
            signs <- intersect(SignatureNames, colnames(tmp))
        } else {
            signs <- Reduce(
                intersect, list(SignatureNames, colnames(tmp), whichSign))}
    } else {stop("There are no signatures in data")}

    tmp <- tmp[,signs]
    tmp <- data.frame(sapply(tmp, range01))

    if(!is.null(groupByAnnot)){
        if(length(groupByAnnot)!=nrow(tmp)){
            stop("groupByAnnot length is different than samples dimension")}
        if(!is.null(selectByAnnot)){
            if(!all(selectByAnnot %in% groupByAnnot)){
                stop("selectByAnnot is not present in groupByAnnot")}}
    } else {if(!is.null(selectByAnnot)){
        stop("selectByAnnot can be used only if groupByAnnot is also provided")}
    }

    if(is.null(whichSign)){n <- ncol(tmp)} else {n <- length(whichSign)}

    if(!is.null(selectByAnnot)){
        tmp <- tmp[groupByAnnot %in% selectByAnnot,]
        groupByAnnot <- groupByAnnot[groupByAnnot %in% selectByAnnot]}

    tmp1 <- do.call(rbind, lapply(seq_len(ncol(tmp)), function(x){
        data.frame(signvalue=tmp[,x], signature=colnames(tmp[x]),
                   row.names = NULL)}))

    g <- ggplot(tmp1, aes(x=signvalue, y=signature))
    if(is.null(groupByAnnot)){
        g <- g + geom_density_ridges(alpha=0.5, bandwidth = 0.05, scale = 1)
    } else {
        g <- g + geom_density_ridges(aes(fill = rep(groupByAnnot, n)),
                                     alpha=0.5, bandwidth = 0.05, scale = 1)}
    g <- g + scale_fill_discrete(name = "Annotation")
    return(g)
}
