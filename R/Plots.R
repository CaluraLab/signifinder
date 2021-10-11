
SignatureNames <- c("Epithelial", "Mesenchymal", "Piroptosis", "Ferroptosis", "LipidMetabolism",
                    "Hypoxia", "PlatinumResistanceUp", "PlatinumResistanceDown", "Prognostic",
                    "ImmunoScore", "IMR_consensus", "DIF_consensus", "PRO_consensus", "MES_consensus",
                    "IPS", "MHC", "CP", "EC", "SC", "Matrisome", "MitoticIndex")

mycol <- c("#FCFDD4", rev(viridis::magma(10)))
mycol1 <- rev(viridis::viridis(10))
my_colors <- RColorBrewer::brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

GetGenes <- function(name){
    if(name %in% c("Epithelial", "Mesenchymal")){
        g <- EMTdata$Gene_Symbol[EMTdata$Category==name]
    } else if (name %in% c("PlatinumResistanceUp", "PlatinumResistanceDown")){
        g <- PlatinumResistancedata[[name]]
    } else if (name %in% c("IMR_consensus", "DIF_consensus", "PRO_consensus", "MES_consensus")){
        stop("Genes for IMR_consensus, DIF_consensus, PRO_consensus and MES_consensus are not available")
    } else if(name %in% c("MHC", "CP", "EC", "SC")){
        g <- IPSdata$GENE[IPSdata$CLASS==name]
    } else {
        datavar <- eval(parse(text = paste0(name, "data")))
        if(name %in% c("Ferroptosis", "Hypoxia", "ImmunoScore", "IPS", "LipidMetabolism", "Piroptosis")){g <- datavar[,1]
        } else if (name %in% c("Prognostic")){g <- names(datavar$Genes)
        } else if (name %in% c("Matrisome", "MitoticIndex")){g <- datavar}
    }
    res <- cbind(g, rep(name, length(g)))
    colnames(res) <- c("Gene", "Signature")
    return(res)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

signatureNameCheck <- function(data, sName){
    if(!all(sName %in% SignatureNames)){
        stop(paste("signatures must be among:", paste(SignatureNames, collapse = ", ")))}
    if(!all(sName %in% colnames(SummarizedExperiment::colData(data)))){
        stop("signature names must be in data")}
}

#' Scatterplot for a single signature
#'
#' Given a signatures it returns a scatterplot.
#'
#' @param data output from a signature function
#' @param whichSign the signature to plot. It has to be only one
#' @param statistics It should be "mean", "median" or "quantiles" to be plot in the graph.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median quantile
#'
#' @export
oneSignPlot <- function(data, whichSign, statistics = NULL){

    if(length(whichSign)>1){stop("you must provide only one signature for this plot")}

    signatureNameCheck(data, whichSign)

    if(!(is.null(statistics))){
        if(!(statistics %in% c("mean", "median", "quantiles"))){
            stop("Statistics to be shown in the graphs must be one of: mean, median and quantiles.")}}

    if(is.vector(data)){signval <- sort(data)
    } else if(is.data.frame(data)){signval <- sort(as.vector(data[whichSign,]))
    } else {signval <- sort(data.frame(colData(data))[,whichSign])}

    g1 <- ggplot() +
        geom_point(mapping = aes(signval, seq_along(signval)), size = 1) +
        labs(x = whichSign, y = "Index") +
        theme(panel.background = element_blank(), axis.line = element_line(colour = "grey50"))
    g2 <- ggplot() +
        geom_histogram(mapping = aes(signval, ..density..), colour = "white", fill = "skyblue2") +
        geom_density(mapping = aes(signval), size = 1.5) +
        labs(x = whichSign, y = "Density") +
        theme(panel.background = element_blank(), axis.line = element_line(colour = "grey50"))

    if(!is.null(statistics)){
        if(statistics=="mean"){
            g1 <- g1 + geom_vline(mapping = aes(xintercept = mean(signval)), col = "red", size = 1.2, linetype = 2)
            g2 <- g2 + geom_vline(mapping = aes(xintercept = mean(signval)), col = "red", size = 1.2, linetype = 2)}
        else if (statistics=="median"){
            g1 <- g1 + geom_vline(mapping = aes(xintercept = median(signval)), col = "red", size = 1.2, linetype = 2)
            g2 <- g2 + geom_vline(mapping = aes(xintercept = median(signval)), col = "red", size = 1.2, linetype = 2)}
        else if (statistics=="quantiles"){
            g1 <- g1 + geom_vline(mapping = aes(xintercept = quantile(signval, 0.25)), col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.50)), col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.75)), col = "red", size = 1.2, linetype = 2)
            g2 <- g2 + geom_vline(mapping = aes(xintercept = quantile(signval, 0.25)), col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.50)), col = "red", size = 1.2, linetype = 2) +
                geom_vline(mapping = aes(xintercept = quantile(signval, 0.75)), col = "red", size = 1.2, linetype = 2)}}
    g1 <- g1  + annotate("text", label = statistics, x = quantile(signval, 0.10), y = length(signval), size = 4, colour = "red")

    return(g1 + g2)
}


#' Scatterplot for a single signature
#'
#' Given a signatures...
#'
#' @param data output from a signature function
#' @param whichSign indicare quali signature usare, di cui verranno presi i geni
#' @param splitBySign se splittare o no le righe dell'heatmap in base alla signature di provenienza dei geni
#' @param sampleAnnot un vettore di annotazione dei sample
#' @param splitByAnnot se splittare o no le colonne in base all'annotazione dei sample, in
#' questo caso sampleAnnot deve essere una variabile categorica
#'
#' @return A ComplexHeatmap object
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation '%v%' HeatmapAnnotation
#' @importFrom magrittr '%>%'
#' @importFrom dplyr group_by summarise_all
#' @importFrom SummarizedExperiment colData
#'
#' @export
geneHeatmapSignPlot <- function(data, whichSign, splitBySign = FALSE,
                                sampleAnnot = NULL, splitByAnnot = FALSE){

    signatureNameCheck(data, whichSign)

    dataset <- getMatrix(data)

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=ncol(dataset)){stop("sampleAnnot length is different than samples dimension")}
    } else {if(splitByAnnot){stop("splitByAnnot can be TRUE only if sampleAnnot is provided")}}

    if(is.vector(data)){
        signval <- matrix(data, nrow = 1, dimnames = list(attr(data, "Signature Name"), colnames(dataset)))
    } else if(is.data.frame(data)){
        signval <- data[whichSign,]
        if(length(whichSign)==1){
            signval <- matrix(signval, nrow = 1, dimnames = list(whichSign, colnames(dataset)))
        } else {
            signval <- sapply(t(signval), range01)
            signval <- as.matrix(t(signval))}
    } else {
        signval <- colData(data)[,whichSign]
        if(length(whichSign)==1){
            signval <- matrix(signval, nrow = 1, dimnames = list(whichSign, colnames(dataset)))
        } else {
            signval <- sapply(signval, range01)
            signval <- as.matrix(t(signval))}
        }

    geneTable <- as.data.frame(do.call(rbind, lapply(whichSign, GetGenes))) %>%
        group_by(Gene) %>% summarise_all(paste, collapse=",")
    signatureGenes <- geneTable$Gene

    filtdataset <- as.matrix(dataset[row.names(dataset) %in% signatureGenes, ])

    htargs <- list(matrix = log2(filtdataset+1), name = "Genes", show_column_names = F, col = mycol,
                   row_names_gp = grid::gpar(fontsize = 5))

    if(length(whichSign)!=1){
        signAnnot <- geneTable$Signature[geneTable$Gene %in% rownames(filtdataset)]
        if(splitBySign){
            htargs$row_split = signAnnot
        } else {
            ha <- rowAnnotation(signature = signAnnot)
            htargs$right_annotation = ha}}

    if(splitByAnnot & is.character(sampleAnnot)){
        htargs$column_split = sampleAnnot
        ht <- Heatmap(signval, name = "Signature", col = mycol1, column_split = sampleAnnot)
    } else {
        if(!is.null(sampleAnnot)){
            hatop = HeatmapAnnotation(sampleAnnot = sampleAnnot)
            htargs$top_annotation = hatop}
        ht <- Heatmap(signval, name = "Signature", col = mycol1)}

    ht2 <- do.call(Heatmap, htargs)
    g <- ht %v% ht2

    return(g)
}


#' Heatmap globale o con una signature che guida il clustering
#'
#' prende in input una matrice di valori che ha le signatures sulle righe e i pazienti sulle colonne.
#'
#' @param data output from a signature function
#' @param whichSign le signature che si vuole far vedere nell'heatmap
#' @param clusterBySign una signature (o più) che guidi il clustering delle colonne dell'heatmap
#' @param sampleAnnot un vettore di annotazione dei sample
#' @param splitByAnnot se splittare o no le colonne in base all'annotazione dei sample, in
#' questo caso sampleAnnot deve essere una variabile categorica
#' @param ... Heatmap params
#'
#' @return A ComplexHeatmap object
#'
#' @importfrom ComplexHeatmap Heatmap '%v%'
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
        if(length(sampleAnnot)!=nrow(data)){stop("sampleAnnot length is different than samples dimension")}
    } else {if(splitByAnnot){stop("splitByAnnot can be TRUE only if sampleAnnot is provided")}}

    if(!is.null(whichSign)){data <- data[, intersect(colnames(data), c(whichSign, clusterBySign))]}

    data <- sapply(data, range01)
    data <- as.matrix(t(data))

    dots <- list(...)
    htargs <- matchArguments(dots, list(name = "Signatures", show_column_names = F, col = mycol))

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
            ht <- Heatmap(fm, name = "Guiding Signatures", col = mycol1, column_split = sampleAnnot)
        } else {
            ht <- Heatmap(fm, name = "Guiding Signatures", col = mycol1)}
        g <- ht %v% do.call(Heatmap, htargs)
    }
    return(g)
}


#' Correlation Plot
#'
#'
#'
#' @param data output from a signature function
#' @param whichSign le signature che si vuole far vedere nel correlation plot
#' @param sampleAnnot un vettore di annotazione dei sample, deve essere categorico
#' @param selectByAnnot a group from sampleAnnot to use to construct the correlation plot
#'
#' @return A correlation ellipse graph
#'
#' @importFrom openair corPlot
#' @importFrom SummarizedExperiment colData
#'
#' @export
correlationSignPlot <- function(data, whichSign = NULL, sampleAnnot = NULL, selectByAnnot = NULL){

    if(!is.null(whichSign)){signatureNameCheck(data, whichSign)}

    if(is.data.frame(data)){tmp <- data} else {tmp <- colData(data)}

    if(sum(colnames(tmp) %in% SignatureNames)>0){
        if(is.null(whichSign)){
            signs <- intersect(SignatureNames, colnames(tmp))
        } else {
            signs <- Reduce(intersect, list(SignatureNames, colnames(tmp), whichSign))}
    } else {stop("There are no signatures in data")}

    tmp <- tmp[,signs]

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=nrow(tmp)){stop("sampleAnnot length is different than samples dimension")}
        if(!is.null(selectByAnnot)){
            if(!(selectByAnnot %in% sampleAnnot)){stop("selectByAnnot is not present in sampleAnnot")}
        } else {stop("sampleAnnot can be used only if selectByAnnot is also provided")}
    } else {if(!is.null(selectByAnnot)){stop("selectByAnnot can be used only if sampleAnnot is also provided")}}

    if(!is.null(sampleAnnot)){if(!is.null(selectByAnnot)){tmp <- tmp[sampleAnnot==selectByAnnot,]}}

    SignMatrix <- sapply(tmp, range01)

    # corsign <- cor(as.matrix(SignMatrix))
    # ord <- order(corsign[1, ])
    # corsign_ord <- corsign[ord, ord]
    # g <- ellipse::plotcorr(corsign_ord, col = my_colors[corsign_ord*50+50], mar = c(1,1,1,1))

    g <- corPlot(as.data.frame(SignMatrix), cluster = T, dendrogram = T, lower = T)

    return(g)
}


#' Survival Plot
#'
#'
#'
#' @param data output from a signature function
#' @param survData dataframe che deve avere due colonne, la prima con i dati di survival (tempo) e la seconda con lo
#' status, inoltre le righe devono essere nominate con i nomi dei sample, prenderlo dalla documentazione della KM
#' @param whichSign la signatura di cui si vuol testare la survival
#' @param cutpoint documentazione originale KM
#' @param sampleAnnot deve essere un vettore nominato con i nomi dei sample, come le righe, deve essere categorico
#' @param selectByAnnot a group from sampleAnnot to use to compute the survival
#'
#' @return al momento ritorna sia il plot che delle statistiche sui dati
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median
#'
#' @export
survivalSignPlot <- function(data, survData, whichSign, cutpoint = "mean",
                             sampleAnnot = NULL, selectByAnnot = NULL){

    if(length(whichSign)>1){stop("you must provide only one signature for this plot")}
    signatureNameCheck(data, whichSign)

    if(ncol(survData)>2){stop("survData should contain only two columns with survival data and status information")}

    if(!(is.numeric(cutpoint))){
        if(!(cutpoint %in% c("mean", "median", "optimal"))){
            stop("Cutpoint must be either a number or one of: mean, median and optimal")}}

    if(is.vector(data)){
        tmp <- intersect(names(data), rownames(survData))
        tmp <- as.data.frame(cbind(data[tmp], survData[tmp,]))
    } else if(is.data.frame(data)){
        tmp <- intersect(rownames(data), rownames(survData))
        tmp <- as.data.frame(cbind(data[tmp,], survData[tmp,]))
    } else {
        tmp <- intersect(rownames(colData(data)), rownames(survData))
        tmp <- as.data.frame(cbind(colData(data)[tmp,], survData[tmp,]))}

    grp  <- rep("high", nrow(tmp))
    names(grp) <- rownames(tmp)
    n <- ncol(tmp)
    colnames(tmp)[c(n-1,n)] <- c("survival", "status")
    if(cutpoint=="mean"){
        grp[which(tmp[,whichSign] < mean(tmp[,whichSign]))] <- "low"
    } else if(cutpoint=="median"){
        grp[which(tmp[,whichSign] < median(tmp[,whichSign]))] <- "low"
    } else if(cutpoint=="optimal"){
        optval <- maxstat::maxstat.test(survival::Surv(survival, status) ~ tmp[,whichSign],
                                        data = tmp, smethod="LogRank", pmethod="Lau94")
        grp[which(tmp[,whichSign] < optval$estimate)] <- "Low"
    } else {grp[which(tmp[,whichSign] < cutpoint)] <- "low"}

    if((sum(grp=="low")<length(grp)/10) | (sum(grp=="low")>length(grp)*9/10)){
        warning(paste("groups size is non homogeneous:", sum(grp=="low"), "low and", sum(grp=="high"), "high"))}

    if(!is.null(sampleAnnot)){
        if(length(sampleAnnot)!=nrow(tmp)){stop("sampleAnnot length is different than samples dimension")}
        if(!is.null(selectByAnnot)){
            if(!(selectByAnnot %in% sampleAnnot)){stop("selectByAnnot is not present in sampleAnnot")}
        } else {stop("sampleAnnot can be used only if selectByAnnot is also provided")}
    } else {if(!is.null(selectByAnnot)){stop("selectByAnnot can be used only if sampleAnnot is also provided")}}

    tmp <- cbind(tmp, grp)
    if(!is.null(sampleAnnot)){if(!is.null(selectByAnnot)){tmp <- tmp[sampleAnnot==selectByAnnot,]}}

    fit <- survival::survfit(survival::Surv(survival, status) ~ grp, data = tmp)

    g <- survminer::ggsurvplot(fit, data = tmp, risk.table = T, legend.title = whichSign,
                               palette = c("red", "blue"), ggtheme = ggplot2::theme_gray(15),
                               font.legend = 15, font.tickslab = 15, font.x = 15, font.y = 15,
                               risk.table.fontsize = 5, pval = T, surv.median.line = "hv",
                               risk.table.col = "strata")
    return(g)
}


#' Ridgeline Plot
#'
#'
#'
#' @param data output from a signature function
#' @param whichSign le signature che si vuole far vedere nel ridgeline plot
#' @param groupByAnnot un vettore di annotazione dei sample, deve essere categorico
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom SummarizedExperiment colData
#'
#' @export
ridgelineSignPlot <- function(data, whichSign = NULL, groupByAnnot = NULL){

    if(!is.null(whichSign)){signatureNameCheck(data, whichSign)}

    if(is.data.frame(data)){tmp <- data} else {tmp <- colData(data)}

    if(sum(colnames(tmp) %in% SignatureNames)>0){
        if(is.null(whichSign)){
            signs <- intersect(SignatureNames, colnames(tmp))
        } else {
            signs <- Reduce(intersect, list(SignatureNames, colnames(tmp), whichSign))}
    } else {stop("There are no signatures in data")}

    tmp <- tmp[,signs]
    tmp <- data.frame(sapply(tmp, range01))

    if(!is.null(groupByAnnot)){if(length(groupByAnnot)!=nrow(tmp)){stop("groupByAnnot length is different than samples dimension")}}

    if(is.null(whichSign)){n <- ncol(tmp)} else {n <- length(whichSign)}

    tmp1 <- do.call(rbind, lapply(seq_len(ncol(tmp)), function(x){
        data.frame(signvalue=tmp[,x], signature=colnames(tmp[x]), row.names = NULL)}))

    g <- ggplot(tmp1, aes(x=signvalue, y=signature))
    if(is.null(groupByAnnot)){
        g <- g + geom_density_ridges(alpha=0.5)
    } else {
        g <- g + geom_density_ridges(aes(fill = rep(groupByAnnot, n)), alpha=0.5)}
    return(g)
}
