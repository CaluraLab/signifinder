SignatureNames <- c(
    "EMT_Miow-Epithelial", "EMT_Miow-Mesenchymal", "EMT_Mak", "EMT_Cheng",
    "Pyroptosis_Ye", "Pyroptosis_Shao", "Pyroptosis_Lin", "Pyroptosis_Li",
    "Ferroptosis_Liang", "Ferroptosis_Li", "Ferroptosis_Liu", "Ferroptosis_Ye",
    "LipidMetabolism_Zheng",
    "Hypoxia_Buffa",
    "PlatinumResistanceUp", "PlatinumResistanceDown",
    "ImmunoScore_Hao", "ImmunoScore_Roh",
    "ConsensusOV_Chen-IMR", "ConsensusOV_Chen-DIF",
    "ConsensusOV_Chen-PRO", "ConsensusOV_Chen-MES",
    "IPS_Charoentong", "IPS_Charoentong-MHC", "IPS_Charoentong-CP",
    "IPS_Charoentong-EC", "IPS_Charoentong-SC",
    "Matrisome_Yuzhalin",
    "MitoticIndex_Yang",
    "ImmuneCyt_Rooney", "ImmuneCyt_Davoli",
    "IFN_Ayers", "ExpandedImmune_Ayers", "Tinflam_Ayers",
    "TLS",
    "StemCellCD49f_Smith",
    "CIN_Carter-25", "CIN_Carter-70",
    "CellCycle_Lundberg", "CellCycle_Davoli",
    "Chemokines_Messina",
    "ASC_Smith",
    "PassON_Du",
    "IPRES",
    "CIS",
    "Glycolysis_Zhang", "Glycolysis_Xu",
    "Autophagy_Xu", "Autophagy_Wang", "Autophagy_ChenM-OS",
    "Autophagy_ChenM-DFS", "Autophagy_ChenH",
    "ECM_Chakravarthy_up", "ECM_Chakravarthy_down",
    "HRDS_Lu",
    "ISC_MerlosSuarez-ISCEphB2", "ISC_MerlosSuarez-LateTA",
    "ISC_MerlosSuarez-ISCLgr5", "ISC_MerlosSuarez-Prolif",
    "VEGF_Hu",
    "DNArep_Kang",
    "IPSOV")

mycol <- c("#FCFDD4", rev(viridis::magma(10)))
mycol1 <- rev(viridis::viridis(10))
my_colors <- RColorBrewer::brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

GetGenes <- function(name){
    if(name %in% c("EMT_Miow-Epithelial", "EMT_Miow-Mesenchymal")){
        name <- paste0(substring(name, 10), "-like")
        g <- EMT_Miow$SYMBOL[EMT_Miow$class==name]
    } else if (name %in% c("PlatinumResistanceUp", "PlatinumResistanceDown")){
        g <- PlatinumResistancedata$Gene_Symbol[
            PlatinumResistancedata$Category==name]
    } else if(name %in% c("ECM_Chakravarthy_up", "ECM_Chakravarthy_down")){
        name <- substring(name, 18)
        g <- ECM_Chakravarthy$SYMBOL[grepl(name, ECM_Chakravarthy$class)]
    } else if(name %in% c("CISup", "CISdown")){
        g <- CISdata$Gene_Symbol[CISdata$Category==name]
    } else if (name %in% c(
        "ConsensusOV_Chen-IMR", "ConsensusOV_Chen-DIF",
        "ConsensusOV_Chen-PRO", "ConsensusOV_Chen-MES")){
        name <- substring(name, 18)
        g <- ConsensusOV_Chen$SYMBOL[ConsensusOV_Chen$class==name]
    } else if(name %in% c("IPS_Charoentong-MHC", "IPS_Charoentong-CP",
                          "IPS_Charoentong-EC", "IPS_Charoentong-SC")){
        g <- IPS_Charoentong$SYMBOL[IPS_Charoentong$class==substring(name, 17)]
    } else if(name %in% c("CIN_Carter-25", "CIN_Carter-70")){
        g <- if(name == "CIN_Carter-25"){
            CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"]
             } else {CIN_Carter$SYMBOL}
    } else if (name %in% c(
        "ISC_MerlosSuarez-ISCEphB2", "ISC_MerlosSuarez-LateTA",
        "ISC_MerlosSuarez-ISCLgr5", "ISC_MerlosSuarez-Prolif")){
        name <- substring(name, 18)
        g <- ISC_MerlosSuarez$SYMBOL[ISC_MerlosSuarez$class == name]
    } else if(name %in% "IPRES"){
        g <- as.vector(unlist(IPRESdata))
    } else if(name == "Tinflam_Ayers"){
        g <- Tinflam_Ayers$SYMBOL[Tinflam_Ayers$class=="TInflam"]
    } else if(name %in% c("Autophagy_ChenM-OS", "Autophagy_ChenM-DFS")){
        name <- substring(name, 17)
        g <- Autophagy_ChenM$SYMBOL[Autophagy_ChenM$class == name]
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
        "CellCycle_Lundberg")){
        g <- unique(eval(parse(text = name))[,"SYMBOL"])
    } else {
        datavar <- eval(parse(text = paste0(name, "data")))
        if(name %in% c("IPSOV")){
            g <- datavar[,1]
        } else if (name %in% c("TLS")){
            g <- datavar}
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
                assays = datasetm, colData = data.frame(name = result))
            colnames(SummarizedExperiment::colData(result)) <- SignName
        } else {
            result <- SummarizedExperiment::SummarizedExperiment(
                assays = datasetm, colData = t(result))
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
        cat("Performing permutation number", x, "\n")
        permlist <- lapply(seq_len(length(gset.idx.list)), function(i)
            sample(datasetGenes, size = lengths(filteredGeneSets)[i],
                   replace = FALSE))
        args$gset.idx.list <- permlist
        gsva_matrix <- suppressWarnings(do.call(gsva, args))
        data.frame(t(gsva_matrix))}, mc.cores = 1)
    permutedResByGeneSet <- split.default(x = Reduce(cbind, permutedResults),
                                          seq_len(length(gset.idx.list)))
    permutedResByGeneSet <- lapply(
        permutedResByGeneSet, function(x)data.frame(t(x)))
    finalRes <- do.call(rbind, lapply(
        seq_len(length(gset.idx.list)), function(i){
            gspvalues <- sapply(1:ncol(expr), function(j){
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

coefficientsScore <- function(ourdata, datasetm, nametype, namesignature){

    ourdata$Gene_Symbol <- geneIDtrans(nametype, ourdata$Gene_Symbol)

    percentageOfGenesUsed(namesignature, datasetm, ourdata$Gene_Symbol)

    ourdata <- ourdata[ourdata$Gene_Symbol %in% row.names(datasetm), ]
    columnNA <- managena(datasetm = datasetm, genes = ourdata$Gene_Symbol)
    score <- colSums(datasetm[ourdata$Gene_Symbol, ] * ourdata$Coefficient,
                     na.rm = TRUE)
    score[columnNA > 0.9] <- NA
    return(score)
}

coeffScore <- function(sdata, datasetm, namesignature, detail = NULL){

    percentageOfGenesUsed(namesignature, datasetm, sdata$SYMBOL, detail)

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
                      namesignature){

    genes <- geneIDtrans(nametype, genes)

    percentageOfGenesUsed(namesignature, datasetm, genes)

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

dataTransformation <- function(data, trans, hgReference){

    if(hgReference=="hg19"){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if(hgReference=="hg38"){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {stop("Human reference genome must be either hg19 or hg38")}

    exons.db = ensembldb::exonsBy(txdb, by = "gene")

    g <- rownames(data)
    egs <- suppressMessages(mapIds(org.Hs.eg.db, keys = g, column = "ENTREZID",
                                   keytype = "SYMBOL", multiVals = "first"))
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

    return(tdata)
}

percentageOfGenesUsed <- function(name, datasetm, gs, detail = NULL){

    g_per <- (sum(gs %in% row.names(datasetm))/length(gs)) * 100

    if(is.null(detail)){
        message(paste0(name, " function is using ", round(g_per), "% of genes"))
        if(g_per < 30){
            warning("The signature is computed with less than 30% of its genes")
        }
    } else {
        message(paste0(name, " function is using ", round(g_per),
                       "% of ", detail, " genes"))
        if(g_per < 30){
            warning(paste(detail,"is computed with less than 30% of its genes"))
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


#' 7 spots resolution
#'
#' Given a 10X Visium dataset, it reassigns to each spot the aggregation of it
#' with the nearest.
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
    spotcoords <- data.frame(
        row = as.vector(sapply(seq(0,76,2), function(i) rep(c(i,i+1),64) )),
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
#' It shows a table containing all the information of the signatures collected
#' in the package.
#'
#' @param tumorTissue type of tissue for which the signature is developed
#' @param signatureType filter for signature type
#' @param inputType type of data you are using: microarray or rnaseq
#' @param description whether to show or not the signature description
#'
#' @return a data frame
#'
#' @export
availableSignatures <- function(tumorTissue = NULL, signatureType = NULL,
                                inputType = NULL, description = TRUE){
    st <- signatureTable
    if(!is.null(tumorTissue)){st <- st[st$tumorTissue==tumorTissue,]}
    if(!is.null(signatureType)){st <- st[st$signatureType==signatureType,]}
    if(!is.null(inputType)){st <- st[st$inputType==inputType,]}
    if(!description){st <- st[,1:6]}
    return(st)
}
