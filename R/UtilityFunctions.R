SignatureNames <- c(
    "EMT_Miow_Epithelial", "EMT_Miow_Mesenchymal", "EMT_Mak", "EMT_Cheng",
    "Pyroptosis_Ye", "Pyroptosis_Shao", "Pyroptosis_Lin", "Pyroptosis_Li",
    "Ferroptosis_Liang", "Ferroptosis_Li", "Ferroptosis_Liu", "Ferroptosis_Ye",
    "LipidMetabolism_Zheng",
    "Hypoxia_Buffa",
    "ImmunoScore_Hao", "ImmunoScore_Roh",
    "ConsensusOV_Chen_IMR", "ConsensusOV_Chen_DIF",
    "ConsensusOV_Chen_PRO", "ConsensusOV_Chen_MES",
    "IPS_Charoentong", "IPS_Charoentong_MHC", "IPS_Charoentong_CP",
    "IPS_Charoentong_EC", "IPS_Charoentong_SC",
    "Matrisome_Yuzhalin",
    "MitoticIndex_Yang",
    "ImmuneCyt_Rooney", "ImmuneCyt_Davoli",
    "IFN_Ayers", "ExpandedImmune_Ayers", "Tinflam_Ayers",
    "TLS_Cabrita",
    "StemCellCD49f_Smith",
    "CIN_Carter_25", "CIN_Carter_70",
    "CellCycle_Lundberg", "CellCycle_Davoli",
    "Chemokines_Messina",
    "ASC_Smith",
    "PassON_Du",
    "IPRES_Hugo",
    "CIS_Robertson",
    "Glycolysis_Zhang", "Glycolysis_Xu",
    "Autophagy_Xu", "Autophagy_Wang", "Autophagy_ChenM_OS",
    "Autophagy_ChenM_DFS", "Autophagy_ChenH",
    "ECM_Chakravarthy_up", "ECM_Chakravarthy_down",
    "HRDS_Lu",
    "ISC_MerlosSuarez_ISCEphB2", "ISC_MerlosSuarez_LateTA",
    "ISC_MerlosSuarez_ISCLgr5", "ISC_MerlosSuarez_Prolif",
    "VEGF_Hu",
    "DNArep_Kang",
    "IPSOV")

mycol <- c("#FCFDD4", rev(viridis::magma(10)))
mycol1 <- rev(viridis::viridis(10))
my_colors <- RColorBrewer::brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

GetGenes <- function(name){
    if(name %in% c("EMT_Miow_Epithelial", "EMT_Miow_Mesenchymal")){
        sname <- paste0(substring(name, 10), "-like")
        g <- EMT_Miow$SYMBOL[EMT_Miow$class==sname]
    } else if(name %in% c("ECM_Chakravarthy_up", "ECM_Chakravarthy_down")){
        sname <- substring(name, 18)
        g <- ECM_Chakravarthy$SYMBOL[grepl(sname, ECM_Chakravarthy$class)]
    } else if (name %in% c(
        "ConsensusOV_Chen_IMR", "ConsensusOV_Chen_DIF",
        "ConsensusOV_Chen_PRO", "ConsensusOV_Chen_MES")){
        sname <- substring(name, 18)
        g <- ConsensusOV_Chen$SYMBOL[ConsensusOV_Chen$class==sname]
    } else if(name %in% c("IPS_Charoentong_MHC", "IPS_Charoentong_CP",
                          "IPS_Charoentong_EC", "IPS_Charoentong_SC")){
        g <- IPS_Charoentong$SYMBOL[IPS_Charoentong$class==substring(name, 17)]
    } else if(name %in% c("CIN_Carter_25", "CIN_Carter_70")){
        g <- if(name == "CIN_Carter_25"){
            CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"]
             } else {CIN_Carter$SYMBOL}
    } else if (name %in% c(
        "ISC_MerlosSuarez_ISCEphB2", "ISC_MerlosSuarez_LateTA",
        "ISC_MerlosSuarez_ISCLgr5", "ISC_MerlosSuarez_Prolif")){
        sname <- substring(name, 18)
        g <- ISC_MerlosSuarez$SYMBOL[ISC_MerlosSuarez$class == sname]
    } else if(name == "Tinflam_Ayers"){
        g <- Tinflam_Ayers$SYMBOL[Tinflam_Ayers$class=="TInflam"]
    } else if(name %in% c("Autophagy_ChenM_OS", "Autophagy_ChenM_DFS")){
        sname <- substring(name, 17)
        g <- Autophagy_ChenM$SYMBOL[Autophagy_ChenM$class == sname]
    } else if(name %in% c(
        "EMT_Mak", "EMT_Cheng",
        "Pyroptosis_Ye", "Pyroptosis_Shao", "Pyroptosis_Lin", "Pyroptosis_Li",
        "Ferroptosis_Liang","Ferroptosis_Li","Ferroptosis_Liu","Ferroptosis_Ye",
        "LipidMetabolism_Zheng", "Hypoxia_Buffa", "Matrisome_Yuzhalin",
        "Chemokines_Messina", "ImmunoScore_Hao", "ImmunoScore_Roh",
        "ImmuneCyt_Rooney", "ImmuneCyt_Davoli", "IFN_Ayers", "HRDS_Lu",
        "ExpandedImmune_Ayers", "CellCycle_Davoli", "PassON_Du", "VEGF_Hu",
        "DNArep_Kang", "ASC_Smith", "IPS_Charoentong", "StemCellCD49f_Smith",
        "Glycolysis_Zhang", "Glycolysis_Xu", "MitoticIndex_Yang",
        "Autophagy_Xu", "Autophagy_Wang", "Autophagy_ChenH",
        "CellCycle_Lundberg", "IPRES_Hugo", "CIS_Robertson", "TLS_Cabrita")){
        g <- unique(eval(parse(text = name))[,"SYMBOL"])
    } else {
        datavar <- eval(parse(text = paste0(name, "data")))
        if(name %in% c("IPSOV")){g <- datavar[,1]}
    }
    res <- cbind(g, rep(name, length(g)))
    colnames(res) <- c("Gene", "Signature")
    return(res)
}

range01 <- function(x){
    (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

signatureNameCheck <- function(data, sName){
    if(!all(sName %in% SignatureNames)){
        stop(paste("signatures must be among:",
                   paste(SignatureNames, collapse = ", ")))}
    if(!all(sName %in% colnames(SummarizedExperiment::colData(data)))){
        stop("signature names must be in data")}
}

matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)}

getMatrix <- function(userdata){
    if(!is.matrix(userdata)) {
        if(is(userdata, "Seurat")){
            if(length(userdata@assays)==1){
                userdata <- as.matrix(userdata@assays$RNA@data)
            } else {userdata <- as.matrix(userdata@assays$SCT@data)}
        } else if(class(userdata) %in% c("SpatialExperiment",
                                         "SummarizedExperiment",
                                         "SingleCellExperiment")){
            userdata <- as.matrix(SummarizedExperiment::assay(userdata))
        } else if(is.data.frame(userdata)){userdata <- as.matrix(userdata)
        } else {stop("This dataset type is not supported")}}
    return(userdata)}

returnAsInput <- function(userdata, result, SignName, datasetm){
    if(!is.matrix(userdata) & !is.data.frame(userdata)) {
        if(is(userdata, "Seurat")){
            if(is.vector(result)){
                names <- c(colnames(userdata@meta.data), SignName)
                userdata@meta.data <- cbind(userdata@meta.data, name=result)
                colnames(userdata@meta.data) <- names
            } else {
                userdata@meta.data <- cbind(userdata@meta.data, t(result))}
        } else if(class(userdata) %in% c("SpatialExperiment",
                                         "SummarizedExperiment",
                                         "SingleCellExperiment")){
            names <- c(colnames(userdata@colData), SignName)
            if(is.vector(result)){
                userdata@colData <- cbind(userdata@colData, name=result)
                colnames(userdata@colData) <- names
            } else {
                userdata@colData <- cbind(userdata@colData, t(result))}}
        return(userdata)
    } else if(is.matrix(userdata) | is.data.frame(userdata)){
        if(is.vector(result)){
            result <- SummarizedExperiment::SummarizedExperiment(
                assays = list(norm_expr=datasetm), colData = data.frame(name = result))
            colnames(SummarizedExperiment::colData(result)) <- SignName
        } else {
            result <- SummarizedExperiment::SummarizedExperiment(
                assays = list(norm_expr=datasetm), colData = t(result))
        }
        return(result)}
}

ipsmap <- function(x){
    if(is.na(x)){NA} else {if(x<=0){0} else if(x>=3){10
    } else {round(x*10/3, digits=0)}}}

GSVAPvalues <- function(expr, gset.idx.list, gsvaResult, nperm, args){
    datasetGenes <- rownames(expr)
    filteredGeneSets <- lapply(gset.idx.list, y = datasetGenes, intersect)
    permutedResults <- parallel::mclapply(seq_len(nperm), function(x){
        message("Performing permutation number", x, "\n")
        permlist <- lapply(seq_len(length(gset.idx.list)), function(i)
            sample(datasetGenes, size = lengths(filteredGeneSets)[i],
                   replace = FALSE))
        args$gset.idx.list <- permlist
        gsva_matrix <- do.call(gsva, args)
        data.frame(t(gsva_matrix))}, mc.cores = 1)
    permutedResByGeneSet <- split.default(x = Reduce(cbind, permutedResults),
                                          seq_len(length(gset.idx.list)))
    permutedResByGeneSet <- lapply(
        permutedResByGeneSet, function(x)data.frame(t(x)))
    finalRes <- do.call(rbind, lapply(
        seq_len(length(gset.idx.list)), function(i){
            gspvalues <- sapply(seq_len(ncol(expr)), function(j){
                (min(c(sum(permutedResByGeneSet[[i]][,j]<=gsvaResult[i,j]),
                       sum(permutedResByGeneSet[[i]][,j]>=gsvaResult[i,j]))
                )+1)/(nperm+1)})
            gspvalues}))
    colnames(finalRes) <- colnames(expr)
    rownames(finalRes) <- paste(names(gset.idx.list), "pval", sep = "_")
    return(finalRes)}

consistencyCheck <- function(nametype, functionName, author = NULL){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL"))){
        stop("The name of genes must be either SYMBOL, ENTREZID or ENSEMBL")}

    if(!is.null(author)){
        if(!(author %in% signatureTable$author[
            signatureTable$functionName==functionName])){
            stop("This author is not present for this signature")}}
}

coeffScore <- function(sdata, datasetm, namesignature,
                       detail = NULL, author = ""){

    percentageOfGenesUsed(namesignature, datasetm, sdata$SYMBOL, detail, author)

    sdata <- sdata[sdata$SYMBOL %in% row.names(datasetm), ]
    columnNA <- managena(datasetm = datasetm, genes = sdata$SYMBOL)
    score <- colSums(datasetm[sdata$SYMBOL, ] * sdata$coeff, na.rm = TRUE)
    score[columnNA > 0.9] <- NA

    return(score)
}

meang <- function(x, na.rm){
    exp(mean(log(x[x>0]), na.rm = na.rm))
}

statScore <- function(genes, datasetm, nametype, typeofstat = "mean",
                      namesignature, author = ""){

    genes <- geneIDtrans(nametype, genes)

    percentageOfGenesUsed(namesignature, datasetm, genes, author = author)

    genes <- genes[genes %in% row.names(datasetm)]

    columnNA <- managena(datasetm = datasetm, genes = genes)
    score <- apply(datasetm[intersect(row.names(datasetm), genes), ],
                   2, typeofstat, na.rm = TRUE)
    score[columnNA > 0.9] <- NA

    return(score)
}

managena <- function(datasetm, genes){
    datasetm <- datasetm[row.names(datasetm) %in% genes, ]
    columnNA <- (length(genes) - colSums(!is.na(datasetm)))/length(genes)
    if(sum(columnNA > 0.9)>0) {
        warning("Some samples in the dataset have more than 90%
                not available (NA) expression values")}
    # datasetm[, columnNA > 0.9] <- NA
    # return(datasetm)
    return(columnNA)
}

dataTransformation <- function(dataset, data, trans, hgReference, nametype){

    if(class(dataset)[1] %in% c(
        "SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
        if(trans %in% names(SummarizedExperiment::assays(dataset))){
            return(dataset) }}

    if(hgReference=="hg19"){
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if(hgReference=="hg38"){
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {stop("Human reference genome must be either hg19 or hg38")}

    exons.db = ensembldb::exonsBy(txdb, by = "gene")

    g <- rownames(data)
    egs <- mapIds(org.Hs.eg.db, keys = g, column = "ENTREZID",
                  keytype = nametype, multiVals = "first")
    data <- data[!is.na(egs),]
    egs <- egs[!is.na(egs)]

    # glen <- getGeneLengthAndGCContent(id = eg, org = "hsa")
    exons_g <- lapply(egs, function(eg){ exons.db[[eg]] })
    data <- data[!sapply(exons_g, is.null),]
    egs <- egs[!sapply(exons_g, is.null)]
    exons_g <- exons_g[!sapply(exons_g, is.null)]
    glen <- sapply(names(egs), function(eg){
        sum(BiocGenerics::width(IRanges::reduce( exons_g[[eg]] ))) })
    tdata <- DGEobj.utils::convertCounts(
        countsMatrix = data, unit = trans, geneLength = glen)

    if(class(dataset)[1] %in% c(
        "SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")) {
        SummarizedExperiment::assays(dataset, withDimnames = FALSE)[[trans]] <- tdata
    } else if(is.matrix(dataset) | is.data.frame(dataset)){
        assays <- list(data, tdata)
        names(assays) <- c("norm_expr", trans)
        dataset <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    }

    return(dataset)
}

percentageOfGenesUsed <- function(name, datasetm, gs, detail = NULL,
                                  author = ""){

    g_per <- (sum(gs %in% row.names(datasetm))/length(gs)) * 100

    if(is.null(detail)){
        message(name, author, " function is using ", round(g_per), "% of genes")
        if(g_per < 30){
            warning("The signature is computed with less than 30% of its genes")
        }
    } else {
        message(name, author, " function is using ", round(g_per),
                "% of ", detail, " genes")
        if(g_per < 30){
            warning(detail,"is computed with less than 30% of its genes")
        }
    }
}

geneIDtrans <- function(nametype, genes){
    if(nametype!="SYMBOL"){
        AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db, keys = genes,
            column = nametype, keytype = "SYMBOL", multiVals = "first")
    } else { genes }
}


#' Show available signatures
#'
#' It shows a table containing all the information of the signatures collected
#' in the package.
#'
#' @param tumor type of tumor for which the signature is developed
#' @param tissue type of tissue for which the signature is developed
#' @param topic filter for signature topic
#' @param requiredInput type of data required: microarray or rnaseq
#' @param description whether to show or not the signature description
#'
#' @return a data frame
#'
#' @export
availableSignatures <- function(tumor = NULL, tissue = NULL,
            topic = NULL, requiredInput = NULL, description = TRUE){
    st <- signatureTable
    if(!is.null(tumor)){
        st <- st[grepl(paste(tumor, collapse = "|"), st$tumor),]}
    if(!is.null(tissue)){
        st <- st[grepl(paste(tissue, collapse = "|"), st$tissue),]}
    if(!is.null(topic)){
        st <- st[grepl(paste(topic, collapse = "|"), st$topic),]}
    if(!is.null(requiredInput)){
        st <- st[grepl(paste(requiredInput, collapse = "|"),st$requiredInput),]}
    if(!description){st <- st[,-10]}
    return(st)
}


#' Multiple Signature Computation
#'
#' @description This function computes all the signatures for a specific
#' 'requiredInput' and for a specific 'tissue' and/or 'topic'.
#'
#' @param dataset Expression values. A data frame or a matrix where rows
#' correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}},
#' \code{\link[SpatialExperiment]{SpatialExperiment}} or
#' \code{\link[SeuratObject]{SeuratObject}} containing an assay
#' where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param inputType type of data you are using: microarray or rnaseq.
#' @param whichSign character vector of signature names to compute.
#' @param tumor character vector of tumor types to compute signatures
#' only from a specific tumor(s) (this can also be pan-cancer).
#' @param tissue character vector of tumor tissue to compute signatures
#' only from a specific tissue (this can also be pan-tissue).
#' @param topic compute signatures only from a specific cancer topic.
#' @param ... other arguments passed on to the signature functions.
#'
#' @return A SummarizedExperiment object in which the signatures' scores
#' are added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @export
multipleSign <- function(dataset, nametype = "SYMBOL",
                         inputType = "rnaseq", whichSign = NULL,
                         tumor = NULL, tissue = NULL,
                         topic = NULL, ...){
    argg <- c(as.list(environment()), list(...))
    st <- availableSignatures(tumor = tumor, tissue = tissue,
        topic = topic, requiredInput = inputType)

    if(!is.null(whichSign)){
        st <- st[grepl(paste(whichSign, collapse = "|"), st$signature),]}

    for(i in seq_len(nrow(st))){
        author <- st[i, "author"]
        argg_i <- c(argg, list(author = author))

        signName <- st[i, "functionName"]
        signFunc <- eval(parse(text = signName))
        signFuncArg <- formals(signFunc)

        argg_i[!(names(argg_i) %in% names(signFuncArg))] <- NULL
        signFuncArg[names(signFuncArg) %in% names(argg_i)] <- NULL
        matchedArg <- c(argg_i, signFuncArg)
        matchedArg[names(matchedArg) == "..."] <- NULL

        dataset <- do.call(signFunc, matchedArg)
        argg$dataset <- dataset
    }
    return(dataset)
}


