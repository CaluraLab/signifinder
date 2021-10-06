matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)}

getMatrix <- function(userdata){
    if(!is.matrix(userdata)) {
        if(class(userdata)=="Seurat"){
            if(length(userdata@assays)==1){userdata <- as.matrix(userdata@assays$RNA@data)
            } else {userdata <- as.matrix(userdata@assays$SCT@data)}
        } else if(class(userdata)%in%c("SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
            userdata <- as.matrix(SummarizedExperiment::assay(userdata))
        } else if(class(userdata)=="data.frame"){userdata <- as.matrix(userdata)
        } else {stop("This dataset type is not supported")}}
    return(userdata)}

returnAsInput <- function(userdata, result, SignName, datasetm){
    if(!is.matrix(userdata) & !is.data.frame(userdata)) {
        if(class(userdata)=="Seurat"){
            if(is.vector(result)){names <- c(colnames(userdata@meta.data), SignName)
            userdata@meta.data <- cbind(userdata@meta.data, name=result)
            colnames(userdata@meta.data) <- names
            } else {userdata@meta.data <- cbind(userdata@meta.data, t(result))}
        } else if(class(userdata)%in%c("SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
            names <- c(colnames(userdata@colData), SignName)
            if(is.vector(result)){userdata@colData <- cbind(userdata@colData, name=result)
            colnames(userdata@colData) <- names
            } else {userdata@colData <- cbind(userdata@colData, t(result))}}
        return(userdata)
    } else if(is.matrix(userdata) | is.data.frame(userdata)){
        if(is.vector(result)){
            result <- SummarizedExperiment::SummarizedExperiment(assays=datasetm, colData=data.frame(name=result))
            colnames(SummarizedExperiment::colData(result)) <- SignName
        } else {
            result <- SummarizedExperiment::SummarizedExperiment(assays = datasetm, colData = t(result))
        }
    return(result)}
}

ipsmap <- function(x){
    if (x<=0) {ips <- 0} else if (x>=3) {ips <- 10} else {ips <- round(x*10/3, digits=0)}
    return(ips)}

GSVAPvalues <- function(expr, gset.idx.list, gsvaResult, nperm, args){
    datasetGenes <- rownames(expr)
    filteredGeneSets <- lapply(gset.idx.list, y = datasetGenes, intersect)
    permutedResults <- parallel::mclapply(seq_len(nperm), function(x){
        cat("Performing permutation number", x, "\n")
        permlist <- lapply(seq_len(length(gset.idx.list)), function(i)
            sample(datasetGenes, size = lengths(filteredGeneSets)[i], replace = F))
        args$gset.idx.list <- permlist
        gsva_matrix <- suppressWarnings(do.call(gsva, args))
        data.frame(t(gsva_matrix))}, mc.cores = 1)
    permutedResByGeneSet <- split.default(x = Reduce(cbind, permutedResults), seq_len(length(gset.idx.list)))
    permutedResByGeneSet <- lapply(permutedResByGeneSet, function(x)data.frame(t(x)))
    finalRes <- do.call(rbind, lapply(seq_len(length(gset.idx.list)), function(i){
        gspvalues <- sapply(1:ncol(expr), function(j){
            (min(c(sum(permutedResByGeneSet[[i]][,j]<=gsvaResult[i,j]),
                   sum(permutedResByGeneSet[[i]][,j]>=gsvaResult[i,j])))+1)/(nperm+1)})
        gspvalues}))
    colnames(finalRes) <- colnames(expr)
    rownames(finalRes) <- paste(names(gset.idx.list), "pval", sep = "_")
    return(finalRes)}

firstCheck <- function(nametype, tumorTissue, functionName){
    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL"))){
        stop("The name of genes must be either SYMBOL, ENTREZID or ENSEMBL")}
    if (!(tumorTissue %in% signatureTable$tumorTissue[signatureTable$functionName==functionName])){
        stop("tumorTissue is not available, check availableSignatures() to see which tissues are included")}
    }


#' 7 spots resolution
#'
#' Given a 10X Visium dataset, it reassigns to each spot the aggregation of it with the nearest.
#'
#' @param dataset Seurat object of a 10X Visium dataset
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom Matrix Matrix
#'
#' @export
GetAggregatedSpot <- function(dataset){
    spotcoords <- data.frame(row = as.vector(sapply(seq(0,76,2), function(i) rep(c(i,i+1),64) )),
                             col = rep(seq(0,127),39))
    overlappingspots <- lapply(1:4992, function(x){
        a <- unlist(spotcoords[x,1]) ##row
        b <- unlist(spotcoords[x,2]) ##col
        data.frame(row = c(a, a-1, a-1, a, a+1, a+1, a),
                   col = c(b, b-1, b+1, b+2, b+1, b-1, b-2))})
    counts <- as.matrix(dataset@assays$SCT@data)
    myrows <- dataset@images$slice1@coordinates$row
    mycols <- dataset@images$slice1@coordinates$col
    for(x in seq_len(length(overlappingspots))){
        ovspots <- overlappingspots[[x]]
        if(sum(myrows==ovspots[1,"row"] & mycols==ovspots[1,"col"])){
            ind <- unlist(sapply(1:7, function(i)
                which(myrows==ovspots[i,"row"] & mycols==ovspots[i,"col"])))
            if(length(ind)!=1){
                kcount <- counts[,ind[1]]
                for(i in 2:length(ind)){kcount <- kcount + counts[,ind[i]]}
                counts[,ind[1]] <- kcount}}}
    spcounts <- Matrix(data = counts, sparse = TRUE)
    dataset@assays$Aggregated
    ## now how to add counts to the seurat object??
    return(dataset)
}


#' Show available signatures
#'
#' It shows a table containing all the information of the signatures collected in the package.
#'
#' @param tissue filter for tissue
#' @param signatureType filter for kind of signature
#' @param description whether to show or not the signature description
#'
#' @return a data frame
#'
#' @export
availableSignatures <- function(tumorTissue, signatureType, description){
    # signatureTable
}

