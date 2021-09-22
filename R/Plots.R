
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

plotFirstCheck <- function(sname){
    if(!all(sname %in% SignatureNames)){
        stop(paste("signatures must be among:", paste(SignatureNames, collapse = ", ")))}
}

#' Scatterplot for a single signature
#'
#' Given a signatures it returns a scatterplot.
#'
#' @param signature If data is NULL, signature must be a numeric vector. When data is provided,
#' signature must be the character string of the column name in data.
#' @param data either a data frame or a SummarizedExperiment.
#' @param statistics It should be "mean", "median" or "quantiles" to be plot on the graph.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
oneSignPlot <- function(data, signatureName, statistics = NULL){

    plotFirstCheck(signatureName)

    if(is.vector(data)){
        if(!(is.null(attr(data, "Signature Name")))){
            if(attr(data, "Signature Name")!=signatureName){
                stop(paste("data and signatureName do not match:",attr(data, "Signature Name"),"&",signatureName))}
        } else {
            stop("You are not providing a signature result vector")}
    } else if (is.data.frame(data)){
        if(!(signatureName%in%rownames(data))){
            stop(paste("data and signatureName do not match:", paste(rownames(data), collapse = ", "),"&",signatureName))}
    } else {
        if(!(signatureName%in%colnames(colData(data)))){
            stop(paste("data and signatureName do not match:", paste(colnames(colData(data)), collapse = ", "),"&",signatureName))}
    }

    if(!(is.null(statistics))){
        if(!(statistics%in%c("mean", "median", "quantiles"))){
            stop("Statistics to be shown in the graphs must be one of: mean, median and quantiles.")}}

    if(is.vector(data)){signval <- sort(data)
    } else if(is.data.frame(data)){signval <- sort(as.vector(data[signatureName,]))
    } else {signval <- sort(data.frame(colData(data))[,signatureName])}

    g1 <- ggplot() +
        geom_point(mapping = aes(signval, seq_along(signval)), size = 1) +
        labs(x = signatureName, y = "Index") +
        theme(panel.background = element_blank(), axis.line = element_line(colour = "grey50"))
    g2 <- ggplot() +
        geom_histogram(mapping = aes(signval, ..density..), colour = "white", fill = "skyblue2") +
        geom_density(mapping = aes(signval), size = 1.5) +
        labs(x = signatureName, y = "Density") +
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
#' @param signature If data is NULL, signature must be a numeric vector. When data is provided,
#' signature must be the character string of the column name in data.
#' @param data either a data frame or a SummarizedExperiment.
#'
#' @return A ComplexHeatmap object
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation '%v%' HeatmapAnnotation
#' @importFrom magrittr '%>%'
#' @importFrom dplyr group_by summarise_all
#'
#' @export
geneHeatmapSignPlot <- function(data, signatureName, dataset, splitGenes = FALSE,
                                group = NULL, splitGroups = FALSE){

    plotFirstCheck(signatureName)

    if(is.vector(data)){
        if(!(is.null(attr(data, "Signature Name")))){
            if(attr(data, "Signature Name")!=signatureName){
                stop(paste("data and signatureName do not match:",attr(data, "Signature Name"),"&",signatureName))}
        } else {
            stop("You are not providing a signature result vector")}
    } else if (is.data.frame(data)){
        if(!all(signatureName%in%rownames(data))){
            stop(paste("data and signatureName do not match:", paste(rownames(data), collapse = ", "),"&",signatureName))}
    } else {
        if(!all(signatureName%in%colnames(colData(data)))){
            stop(paste("data and signatureName do not match:", paste(colnames(colData(data)), collapse = ", "),"&",signatureName))}
    }

    if(is.vector(data)){
        signval <- matrix(data, nrow = 1, dimnames = list(attr(data, "Signature Name"), colnames(dataset)))
    } else if(is.data.frame(data)){
        signval <- data[signatureName,]
        if(length(signatureName)==1){
            signval <- matrix(signval, nrow = 1, dimnames = list(signatureName, colnames(dataset)))
        } else {
            signval <- sapply(t(signval), range01)
            signval <- as.matrix(t(signval))}
    } else {
        signval <- colData(data)[,signatureName]
        if(length(signatureName)==1){
            signval <- matrix(signval, nrow = 1, dimnames = list(signatureName, colnames(dataset)))
        } else {
            signval <- sapply(signval, range01)
            signval <- as.matrix(t(signval))}
        }

    geneTable <- as.data.frame(do.call(rbind, lapply(signatureName, GetGenes))) %>%
        group_by(Gene) %>% summarise_all(paste, collapse=",")
    signatureGenes <- geneTable$Gene

    filtdataset <- as.matrix(dataset[row.names(dataset) %in% signatureGenes, ])

    htargs <- list(matrix = log2(filtdataset+1), name = "Genes", show_column_names = F, col = mycol,
                   row_names_gp = grid::gpar(fontsize = 5))

    if(length(signatureName)!=1){
        signAnnot <- geneTable$Signature[geneTable$Gene %in% rownames(filtdataset)]
        if(splitGenes){
            htargs$row_split = signAnnot
        } else {
            ha <- rowAnnotation(signature = signAnnot)
            htargs$right_annotation = ha}}

    if(splitGroups){
            htargs$column_split = group
            ht <- Heatmap(signval, name = "Signature", col = mycol1, column_split = group)
        } else {
            if(!is.null(group)){
                hatop = HeatmapAnnotation(group = group)
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
#' @param data dataframe con le signature sulle colonne oppure oggetto con i colData
#' @param signatureName una signature (o più) che guidi il clustering
#'
#' @return A ComplexHeatmap object
#'
#' @importfrom ComplexHeatmap Heatmap '%v%'
#'
#' @export
heatmapSignPlot <- function(data, signatureName = NULL, group = NULL, splitGroups = FALSE){

    plotFirstCheck(signatureName)

    if(is.data.frame(data)){
        if(sum(colnames(data)%in%SignatureNames)>0){
            data <- data[, colnames(data)%in%SignatureNames]
        } else {stop("There are no signatures in data")}
    } else {
        if(sum(colnames(colData(data))%in%SignatureNames)>0){
            data <- colData(data)[, colnames(colData(data))%in%SignatureNames]
        } else {stop("There are no signatures in data")}}

    if(!all(signatureName %in% colnames(data))){
        stop("Signature(s) in signatureName must be in the data.")}

    data <- sapply(data, range01)
    data <- as.matrix(t(data))

    htargs <- list(name = "Signatures", show_column_names = F, col = mycol)

    if(!is.null(group)){
        if(splitGroups){
            htargs$column_split = group
        } else {
            hatop = HeatmapAnnotation(group = group)
            htargs$top_annotation = hatop}}

    if(is.null(signatureName)){
        htargs$matrix = data
        g <- do.call(Heatmap, htargs)
    } else {
        n <- which(rownames(data) %in% signatureName)
        fm <- as.matrix(data.frame(data)[n,])
        sm <- as.matrix(data.frame(data)[-n,])
        htargs$matrix = sm
        if(splitGroups){
            ht <- Heatmap(fm, name = "Guiding Signatures", col = mycol1, column_split = group)
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
#' @param data dataframe con le signature sulle colonne oppure oggetto con i colData
#' @param signatureName ..
#' @param group ..
#' @param groupToUse ..
#'
#' @return A correlation ellipse graph
#'
#' @importFrom openair corPlot
#'
#' @export
correlationSignPlot <- function(data, signatureName = NULL, group = NULL, groupToUse = NULL){

    if(!is.null(signatureName)){plotFirstCheck(signatureName)}

    if(is.data.frame(data)){tmp <- data} else {tmp <- colData(data)}

    if(is.null(signatureName)){
        signs <- intersect(SignatureNames, colnames(tmp))
    } else {
        signs <- Reduce(intersect, list(SignatureNames, colnames(tmp), signatureName))}

    tmp <- tmp[,signs]

    if(!is.null(group)){if(!is.null(groupToUse)){tmp <- tmp[group==groupToUse,]}}

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
#' @param data matrice di valori delle signatures
#'
#' @return
#'
#'
#' @export
survivalSignPlot <- function(data, survData, signatureName, cutpoint = "mean", group = NULL, groupToUse = NULL){

    plotFirstCheck(signatureName)

    if(is.vector(data)){
        if(!(is.null(attr(data, "Signature Name")))){
            if(attr(data, "Signature Name")!=signatureName){
                stop(paste("data and signatureName do not match:", attr(data, "Signature Name"),"&",signatureName))}
        } else {
            stop("You are not providing a signature result vector")}
    } else if (is.data.frame(data)){
        if(!(signatureName%in%rownames(data))){
            stop(paste("data and signatureName do not match:", paste(rownames(data), collapse = ", "),"&",signatureName))}
    } else {
        if(!(signatureName%in%colnames(colData(data)))){
            stop(paste("data and signatureName do not match:", paste(colnames(colData(data)), collapse = ", "),"&",signatureName))}
    }

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
    if(cutpoint=="mean"){
        grp[which(tmp[,signatureName] < mean(tmp[,signatureName]))] <- "low"
    } else if(cutpoint=="median"){
        grp[which(tmp[,signatureName] < median(tmp[,signatureName]))] <- "low"
    } else if(cutpoint=="optimal"){
        optval <- maxstat::maxstat.test(survival::Surv(os, status) ~ tmp[,signatureName],
                                        data = tmp, smethod="LogRank", pmethod="Lau94")
        grp[which(tmp[,signatureName] < optval$estimate)] <- "Low"
    } else {grp[which(tmp[,signatureName] < cutpoint)] <- "low"}

    if((sum(grp=="low")<length(grp)/10) | (sum(grp=="low")>length(grp)*9/10)){
        warning(paste("groups size is non homogeneous:", sum(grp=="low"), "low and", sum(grp=="high"), "high"))}

    tmp <- cbind(tmp, grp)
    if(!is.null(group)){if(!is.null(groupToUse)){tmp <- tmp[group==groupToUse,]}}

    fit <- survival::survfit(survival::Surv(os, status) ~ grp, data = tmp)

    g <- survminer::ggsurvplot(fit, data = tmp, risk.table = T, legend.title = signatureName,
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
#' @param data matrice di valori delle signatures
#'
#' @return
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#'
#' @export
ridgelineSignPlot <- function(data, signatureName = NULL, group = NULL){

    if(!is.null(signatureName)){plotFirstCheck(signatureName)}

    if(is.data.frame(data)){
        if(!all(signatureName %in% colnames(data))){
            stop(paste("data and signatureName do not match:", paste(rownames(data), collapse = ", "),"&",signatureName))}
    } else {
        if(!all(signatureName %in% colnames(colData(data)))){
            stop(paste("data and signatureName do not match:", paste(colnames(colData(data)), collapse = ", "),"&",signatureName))}
    }

    if(is.data.frame(data)){tmp <- data} else {tmp <- colData(data)}

    if(is.null(signatureName)){
        signs <- intersect(SignatureNames, colnames(tmp))
    } else {
        signs <- Reduce(intersect, list(SignatureNames, colnames(tmp), signatureName))}

    tmp <- tmp[,signs]
    tmp <- data.frame(sapply(tmp, range01))

    if(is.null(signatureName)){n <- ncol(tmp)} else {n <- length(signatureName)}

    tmp1 <- do.call(rbind, lapply(seq_len(ncol(tmp)), function(x){
        data.frame(signvalue=tmp[,x], signature=colnames(tmp[x]), row.names = NULL)}))

    g <- ggplot(tmp1, aes(x=signvalue, y=signature))
    if(is.null(group)){
        g <- g + geom_density_ridges(alpha=0.5)
    } else {
        g <- g + geom_density_ridges(aes(fill = rep(group, n)), alpha=0.5)}
    return(g)
}
