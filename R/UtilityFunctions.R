SignatureNames <- c(
    "EMT_Miow_Epithelial",
    "EMT_Miow_Mesenchymal",
    "EMT_Mak",
    "EMT_Cheng",
    "EMT_Thompson",
    "EMT_Barkley_cEMT",
    "EMT_Barkley_pEMT",
    "Pyroptosis_Ye",
    "Pyroptosis_Shao",
    "Pyroptosis_Lin",
    "Pyroptosis_Li",
    "Ferroptosis_Liang",
    "Ferroptosis_Li",
    "Ferroptosis_Liu",
    "Ferroptosis_Ye",
    "LipidMetabolism_Zheng",
    "Hypoxia_Buffa",
    "Hypoxia_Barkley",
    "ImmunoScore_Hao",
    "ImmunoScore_Roh",
    "ConsensusOV_Chen_IMR",
    "ConsensusOV_Chen_DIF",
    "ConsensusOV_Chen_PRO",
    "ConsensusOV_Chen_MES",
    "IPS_Charoentong",
    "IPS_Charoentong_MHC",
    "IPS_Charoentong_CP",
    "IPS_Charoentong_EC",
    "IPS_Charoentong_SC",
    "Matrisome_Yuzhalin",
    "MitoticIndex_Yang",
    "ImmuneCyt_Rooney",
    "ImmuneCyt_Davoli",
    "IFN_Ayers",
    "ExpandedImmune_Ayers",
    "Tinflam_Ayers",
    "Tinflam_Thompson",
    "TLS_Cabrita",
    "StemCellCD49f_Smith",
    "CIN_Carter_25",
    "CIN_Carter_70",
    "CellCycle_Lundberg",
    "CellCycle_Davoli",
    "CellCycle_Barkley",
    "Chemokines_Messina",
    "ASC_Smith",
    "PassON_Du",
    "IPRES_Hugo",
    "CIS_Robertson",
    "Glycolysis_Zhang",
    "Glycolysis_Xu",
    "Autophagy_Xu",
    "Autophagy_Wang",
    "Autophagy_ChenM_OS",
    "Autophagy_ChenM_DFS",
    "Autophagy_ChenH",
    "ECM_Chakravarthy_up",
    "ECM_Chakravarthy_down",
    "HRDS_Lu",
    "ISC_MerlosSuarez_ISCEphB2",
    "ISC_MerlosSuarez_LateTA",
    "ISC_MerlosSuarez_ISCLgr5",
    "ISC_MerlosSuarez_Prolif",
    "VEGF_Hu",
    "DNArep_Kang",
    "IPSOV_Shen",
    "State_Neftel_MES2",
    "State_Neftel_MES1",
    "State_Neftel_AC",
    "State_Neftel_OPC",
    "State_Neftel_NPC1",
    "State_Neftel_NPC2",
    "State_Barkley_Alveolar",
    "State_Barkley_Basal",
    "State_Barkley_Squamous",
    "State_Barkley_Glandular",
    "State_Barkley_Ciliated",
    "State_Barkley_AC", 
    "State_Barkley_OPC", 
    "State_Barkley_NPC",
    "State_Tirosh_MITF",
    "State_Tirosh_AXL",
    "Combined_Thompson",
    "APM_Thompson",
    "APM_Wang",
    "MPS_PerezGuijarro",
    "IRG_Yang",
    "TGFB_Mariathasan",
    "ADO_Sidders",
    "MITFlowPTENneg_Cabrita",
    "LRRC15CAF_Dominguez",
    "ICBResponse_Chen_responder",
    "ICBResponse_Chen_nonresponder",
    "COXIS_Bonavita",
    "SCSubtype_Wu_Basal",
    "SCSubtype_Wu_Her2E",
    "SCSubtype_Wu_LumA",
    "SCSubtype_Wu_LumB",
    "Stress_Barkley",
    "Interferon_Barkley",
    "Oxphos_Barkley",
    "Metal_Barkley"
)

#' @importFrom viridis magma viridis
mycol <- c("#FCFDD4", rev(magma(10)))
mycol1 <- rev(viridis(10))

#' @importFrom RColorBrewer brewer.pal
my_colors <- brewer.pal(5, "Spectral")

my_colors <- colorRampPalette(my_colors)(100)

.GetGenes <- function(name) {
    if (name %in% c("EMT_Miow_Epithelial", "EMT_Miow_Mesenchymal")) {
        sname <- paste0(substring(name, 10), "-like")
        g <- EMT_Miow$SYMBOL[EMT_Miow$class == sname]
    } else if (name %in% c("EMT_Barkley_cEMT", "EMT_Barkley_pEMT")) {
        sname <- substring(name, 13)
        g <- PanState_Barkley$SYMBOL[PanState_Barkley$class == sname]
    } else if (name %in% c("ECM_Chakravarthy_up", "ECM_Chakravarthy_down")) {
      sname <- substring(name, 18)
      g <- ECM_Chakravarthy$SYMBOL[grepl(sname, ECM_Chakravarthy$class)]
    } else if (name %in% c(
        "ConsensusOV_Chen_IMR", "ConsensusOV_Chen_DIF",
        "ConsensusOV_Chen_PRO", "ConsensusOV_Chen_MES"
    )) {
        sname <- substring(name, 18)
        g <- ConsensusOV_Chen$SYMBOL[ConsensusOV_Chen$class == sname]
    } else if (name %in% c(
        "State_Neftel_MES2", "State_Neftel_MES1",
        "State_Neftel_AC", "State_Neftel_OPC",
        "State_Neftel_NPC1", "State_Neftel_NPC2"
    )) {
        sname <- substring(name, 14)
        g <- State_Neftel$SYMBOL[State_Neftel$class == sname]
    } else if (name %in% c(
      "State_Barkley_Alveolar", "State_Barkley_Basal", "State_Barkley_Squamous",
      "State_Barkley_Glandular", "State_Barkley_Ciliated", "State_Barkley_AC",
      "State_Barkley_OPC", "State_Barkley_NPC"
    )) {
      sname <- substring(name, 15)
      g <- PanState_Barkley$SYMBOL[PanState_Barkley$class == sname]
    } else if (name %in% c(
      "State_Tirosh_MITF", "State_Tirosh_AXL"
    )) {
      sname <- substring(name, 14)
      g <- State_Tirosh$SYMBOL[State_Tirosh$class == sname]
    } else if (name %in% c(
        "IPS_Charoentong_MHC", "IPS_Charoentong_CP",
        "IPS_Charoentong_EC", "IPS_Charoentong_SC"
    )) {
        g <-IPS_Charoentong$SYMBOL[IPS_Charoentong$class == substring(name, 17)]
    } else if (name %in% c("CIN_Carter_25", "CIN_Carter_70")) {
        g <- if (name == "CIN_Carter_25") {
            CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"]
        } else { CIN_Carter$SYMBOL }
    } else if (name %in% c(
        "ISC_MerlosSuarez_ISCEphB2", "ISC_MerlosSuarez_LateTA",
        "ISC_MerlosSuarez_ISCLgr5", "ISC_MerlosSuarez_Prolif"
    )) {
        sname <- substring(name, 18)
        g <- ISC_MerlosSuarez$SYMBOL[ISC_MerlosSuarez$class == sname]
    } else if (name == "Tinflam_Ayers") {
        g <- Tinflam_Ayers$SYMBOL[Tinflam_Ayers$class == "TInflam"]
    } else if (name %in% c("Autophagy_ChenM_OS", "Autophagy_ChenM_DFS")) {
        sname <- substring(name, 17)
        g <- Autophagy_ChenM$SYMBOL[Autophagy_ChenM$class == sname]
    } else if (name == "Combined_Thompson"){
        name_emt <- .GetGenes("EMT_Thompson")[ , "Gene"]
        name_inf <- .GetGenes("Tinflam_Thompson")[ , "Gene"]
        g <- union(name_emt, name_inf)
    } else if (name %in% c(
      "SCSubtype_Wu_Basal", "SCSubtype_Wu_Her2E",
      "SCSubtype_Wu_LumA", "SCSubtype_Wu_LumB"
    )) {
      sname <- substring(name, 14)
      g <- SCSubtype_Wu$SYMBOL[SCSubtype_Wu$class == sname]
    } else if (name %in% c("ICBResponse_Chen_responder",
                           "ICBResponse_Chen_nonresponder"
    )) {
      sname <- substring(name, 18)
      g <- ICBResponse_Chen$SYMBOL[ICBResponse_Chen$class == sname]
    } else if (name %in% c(
      "Hypoxia_Barkley", "Stress_Barkley", "Interferon_Barkley", 
      "Oxphos_Barkley", "Metal_Barkley")) {
      sname <- substring(name, 1, regexpr("_", name) - 1)
      g <- PanState_Barkley$SYMBOL[PanState_Barkley$class == sname]
    } else if (name == "CellCycle_Barkley") {
      sname <- substring(name, 5, regexpr("_", name) - 1)
      g <- PanState_Barkley$SYMBOL[PanState_Barkley$class == sname]
    } else if (name %in% c(
        "EMT_Mak", "EMT_Cheng", "EMT_Thompson", "Pyroptosis_Ye",
        "Pyroptosis_Shao", "Pyroptosis_Lin", "Pyroptosis_Li",
        "Ferroptosis_Liang", "Ferroptosis_Li", "Ferroptosis_Liu",
        "Ferroptosis_Ye", "LipidMetabolism_Zheng", "Hypoxia_Buffa",
        "Matrisome_Yuzhalin", "Chemokines_Messina", "ImmunoScore_Hao",
        "ImmunoScore_Roh", "ImmuneCyt_Rooney", "ImmuneCyt_Davoli", "IFN_Ayers",
        "HRDS_Lu", "ExpandedImmune_Ayers", "CellCycle_Davoli", "PassON_Du",
        "VEGF_Hu", "DNArep_Kang", "ASC_Smith", "IPS_Charoentong",
        "StemCellCD49f_Smith", "Glycolysis_Zhang", "Glycolysis_Xu",
        "MitoticIndex_Yang", "Autophagy_Xu", "Autophagy_Wang",
        "Autophagy_ChenH", "CellCycle_Lundberg", "IPRES_Hugo", "CIS_Robertson",
        "TLS_Cabrita", "IPSOV_Shen", "Tinflam_Thompson", "APM_Thompson",
        "APM_Wang", "MPS_PerezGuijarro", "IRG_Yang", "TGFB_Mariathasan",
        "ADO_Sidders", "MITFlowPTENneg_Cabrita", "LRRC15CAF_Dominguez",
        "COXIS_Bonavita"
    )) {
        g <- unique(eval(parse(text = name))[, "SYMBOL"])
    }
    res <- cbind(g, rep(name, length(g)))
    colnames(res) <- c("Gene", "Signature")
    return(res)
}

#' @importFrom stats var
.range01 <- function(x) {
    if(!var(x, na.rm = TRUE)){
        y <- rep(0.5, length(x))
        y[is.na(y)] <- NA
        y[is.nan(y)] <- NaN
        y
    } else {
        (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
    }
}

#' @importFrom SummarizedExperiment colData
.signatureNameCheck <- function(data, sName) {
    if (!any(sName %in% SignatureNames)) {
        stop("There are no signatures computed with signifinder in whichSign")
    }
    if (!all(sName %in% colnames(colData(data)))) {
        stop("all signature names in whichSign must be in the colData of data")
    }
}

.matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @import SpatialExperiment
setMethod(".getMatrix", "matrix", function(userdata, ...){return(userdata)})
setMethod(".getMatrix", "data.frame", function(userdata, ...){
    userdata <- as.matrix(userdata)
    return(userdata)})
setMethod(".getMatrix", "SummarizedExperiment", function(userdata, whichAssay){
    userdata <- as.matrix(assay(userdata, whichAssay))
    return(userdata)})
setMethod(".getMatrix", "SingleCellExperiment", function(userdata, whichAssay){
    userdata <- as.matrix(assay(userdata, whichAssay))
    return(userdata)})
setMethod(".getMatrix", "SpatialExperiment", function(userdata, whichAssay){
    userdata <- as.matrix(assay(userdata, whichAssay))
    return(userdata)})

#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<-
setMethod(".returnAsInput", signature(userdata = "matrix", result = "vector"),
    function(userdata, result, SignName, datasetm){
        result <- SummarizedExperiment(
            assays = list(norm_expr = datasetm),
            colData = data.frame(name = result))
        colnames(colData(result)) <- SignName
        return(result) })
setMethod(".returnAsInput", signature(userdata ="data.frame", result ="vector"),
          function(userdata, result, SignName, datasetm){
              result <- SummarizedExperiment(
                  assays = list(norm_expr = datasetm),
                  colData = data.frame(name = result))
              colnames(colData(result)) <- SignName
              return(result) })
setMethod(".returnAsInput", signature(userdata ="matrix", result ="matrix"),
    function(userdata, result, SignName, datasetm){
        result <- SummarizedExperiment(
            assays = list(norm_expr = datasetm),
            colData = t(result))
        return(result) })
setMethod(".returnAsInput",signature(userdata="data.frame",result="matrix"),
          function(userdata, result, SignName, datasetm){
              result <- SummarizedExperiment(
                  assays = list(norm_expr = datasetm),
                  colData = t(result))
              return(result) })
setMethod(".returnAsInput",
    signature(userdata = "SummarizedExperiment", result = "vector"),
    function(userdata, result, SignName, datasetm){
        userdata@colData[, SignName] <- result
        return(userdata) })
setMethod(".returnAsInput",
          signature(userdata = "SingleCellExperiment", result = "vector"),
          function(userdata, result, SignName, datasetm){
              userdata@colData[, SignName] <- result
              return(userdata) })
setMethod(".returnAsInput",
          signature(userdata = "SpatialExperiment", result = "vector"),
          function(userdata, result, SignName, datasetm){
              userdata@colData[, SignName] <- result
              return(userdata) })
setMethod(".returnAsInput",
    signature(userdata = "SummarizedExperiment", result = "matrix"),
    function(userdata, result, SignName, datasetm){
        userdata@colData[, colnames(t(result))] <- t(result)
        return(userdata) })
setMethod(".returnAsInput",
          signature(userdata = "SingleCellExperiment", result = "matrix"),
          function(userdata, result, SignName, datasetm){
              userdata@colData[, colnames(t(result))] <- t(result)
              return(userdata) })
setMethod(".returnAsInput",
          signature(userdata = "SpatialExperiment", result = "matrix"),
          function(userdata, result, SignName, datasetm){
              userdata@colData[, colnames(t(result))] <- t(result)
              return(userdata) })

.ipsmap <- function(x) {
    if (is.na(x)) { NA
    } else {
        if (x <= 0) { 0
        } else if (x >= 3) { 10
        } else { round(x*10/3, digits = 0) }}
}

.consistencyCheck <- function(nametype, functionName, author = NULL) {
    if (!(nametype %in% c("SYMBOL", "ENTREZID", "ENSEMBL"))) {
        stop("The name of genes must be either SYMBOL, ENTREZID or ENSEMBL")}
    if (!is.null(author)) {
        if (!(author %in% signatureTable$author[
            signatureTable$functionName==functionName])) {
            stop("This author is not present for this signature")}}
}

.isMalignantCheck <- function(isMalignant, dataset) {
  if(is.null(isMalignant)){
    stop("isMalignant param is missing")
  } else {
    if(length(isMalignant)!=ncol(dataset)){
      stop("length of isMalignant must be equal to ncol(dataset)")}
    if(!is.logical(isMalignant)){
      stop("isMalignant must be a logical vector")}}
}

.coeffScore <- function(
        sdata, datasetm, namesignature, detail = NULL, author = "") {
    .percentageOfGenesUsed(
        namesignature, datasetm, sdata$SYMBOL, detail, author)
    sdata <- sdata[sdata$SYMBOL %in% row.names(datasetm), ]
    columnNA <- .managena(datasetm = datasetm, genes = sdata$SYMBOL)
    score <- colSums(
      datasetm[sdata$SYMBOL, ,drop = FALSE] * sdata$coeff, na.rm = TRUE)
    score[columnNA > 0.9] <- NA
    return(score)
}

.meang <- function(x, na.rm) { exp(mean(log(x[x > 0]), na.rm = na.rm)) }

.statScore <- function(
        genes, datasetm, nametype, typeofstat = "mean",
        namesignature, author = "") {
    genes <- .geneIDtrans(nametype, genes)
    .percentageOfGenesUsed(namesignature, datasetm, genes, author = author)
    genes <- intersect(row.names(datasetm), genes)

    columnNA <- .managena(datasetm = datasetm, genes = genes)
    score <- apply(datasetm[genes,,drop = FALSE], 2, typeofstat, na.rm = TRUE)
    score[columnNA > 0.9] <- NA
    return(score)
}

.managena <- function(datasetm, genes) {
    datasetm <- datasetm[row.names(datasetm) %in% genes, ]
    columnNA <- if(length(genes)==1){
        is.na(datasetm)
    } else {
        (length(genes) - colSums(!is.na(datasetm))) / length(genes)
    }
    if (sum(columnNA > 0.9) > 0) {
        warning("Some samples have more than 90% NA values")}
    return(columnNA)
}

#' @importFrom SummarizedExperiment assays SummarizedExperiment assays<-
#' @importFrom BiocGenerics width
#' @importFrom IRanges reduce
#' @importFrom DGEobj.utils convertCounts
#' @importFrom AnnotationDbi mapIds
#' @importFrom ensembldb exonsBy
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
.dataTransformation <- function(dataset, data, trans, hgReference, nametype) {
    if (is(dataset)[1] %in% c(
        "SpatialExperiment",
        "SummarizedExperiment",
        "SingleCellExperiment"
    )) { if (trans %in% names(assays(dataset))) { return(dataset) } }

    if (hgReference == "hg19") {
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (hgReference == "hg38") {
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    } else { stop("Human reference genome must be either hg19 or hg38") }

    exons.db <- exonsBy(txdb, by = "gene")

    g <- rownames(data)
    egs <- mapIds(
        org.Hs.eg.db, keys = g, column = "ENTREZID",
        keytype = nametype, multiVals = "first")
    use <- !is.na(egs)

    exons_g <- NULL
    exons_g[use] <- lapply(egs[use], function(eg) { exons.db[[eg]] })
    names(exons_g) <- names(egs)
    use <- !vapply(exons_g, is.null, logical(1))

    glen <- NULL
    glen[use] <- vapply(names(egs[use]), function(eg) {
        sum(width(reduce(exons_g[[eg]])))
    }, numeric(1))

    usec <- colSums(data[use,])>0

    tdata <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    tdata[use,usec] <- convertCounts(
        countsMatrix = data[use,usec],
        unit = trans, geneLength = glen[use])

    if (is(dataset)[1] %in% c(
        "SpatialExperiment",
        "SummarizedExperiment",
        "SingleCellExperiment"
    )) { assays(dataset, withDimnames = FALSE)[[trans]] <- tdata
    } else if (is.matrix(dataset) | is.data.frame(dataset)) {
        assays <- list(data, tdata)
        names(assays) <- c("norm_expr", trans)
        dataset <- SummarizedExperiment(assays = assays)}

    return(dataset)
}

.percentageOfGenesUsed <- function(
        name, datasetm, gs, detail = NULL, author = "") {
    g_per <- (sum(gs %in% row.names(datasetm)) / length(gs)) * 100

    if (is.null(detail)) { message(
        name, author, " is using ", round(g_per), "% of signature genes")
        if (g_per < 30) { warning(
            name, author, " is computed with less than 30% of its genes")}
    } else { message(
        name, author, " is using ", round(g_per),
        "% of ", detail, " signature genes")
        if (g_per < 30) { warning(
            name, author, detail,
            " is computed with less than 30% of its genes")}}
    # return(g_per==0)
}

#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
.geneIDtrans <- function(nametype, genes) {
    if (nametype != "SYMBOL") { mapIds(
        org.Hs.eg.db, keys = genes, column = nametype,
        keytype = "SYMBOL", multiVals = "first")
    } else { genes }
}

#' @importFrom stats rpois rnorm
.fakeData <- function(sname, input = "rnaseq") {
    g <- .GetGenes(sname)
    g <- g[, 1]
    n <- length(g) * 5
    if (input == "rnaseq") {
        rmatrix <- matrix(rpois(n, 100), ncol = 5)
    } else if (input == "microarray") {
        rmatrix <- matrix(rnorm(n, 100, 50), ncol = 5)}
    rownames(rmatrix) <- g
    return(rmatrix)
}


#' Show Available Signatures
#'
#' It returns a table with all the information of the signatures collected
#' in signifinder.
#'
#' @param tumor character vector saying the type of tumors for which signatures
#' are developed. Used to filter the signatures in the table.
#' @param tissue character vector saying the type of tissues for which
#' signatures are developed. Used to filter the signatures in the table.
#' @param topic character vector saying the signature topics. Used to filter
#' the signatures in the table.
#' @param requiredInput character string saying the type of data required in
#' input by the signature. Either one of "microarray", "rnaseq" or "sc". Used to
#' filter the signatures in the table.
#' @param description logical. If TRUE it shows the signature's description.
#'
#' @return A data frame with 12 variables:
#' \describe{
#'   \item{signature}{name of the signature}
#'   \item{scoreLabel}{label of the signature when added inside colData section}
#'   \item{functionName}{name of the function to use to compute the signature}
#'   \item{topic}{main cancer topic of the signature}
#'   \item{tumor}{tumor type for which the signature was developed}
#'   \item{tissue}{tumor tissue for which the signature was developed}
#'   \item{cellType}{cell type for which the signature was developed}
#'   \item{requiredInput}{type of data with which the signature was developed}
#'   \item{transformationStep}{data transformation step performed inside the
#'   function starting from the user's 'normArray' or 'normCounts' data}
#'   \item{author}{first author of the work in which the signature is described}
#'   \item{reference}{reference of the work}
#'   \item{description}{signature description and how to evaluate its score}
#'   ...
#' }
#'
#' @examples
#' availableSignatures()
#'
#' @export
availableSignatures <- function(
        tumor = NULL, tissue = NULL, topic = NULL,
        requiredInput = NULL, description = TRUE) {
    st <- signatureTable
    if (!is.null(tumor)) {
        st <- st[grepl(paste(tumor, collapse = "|"), st$tumor), ]
    }
    if (!is.null(tissue)) {
        st <- st[grepl(paste(tissue, collapse = "|"), st$tissue), ]
    }
    if (!is.null(topic)) {
        st <- st[grepl(paste(topic, collapse = "|"), st$topic), ]
    }
    if (!is.null(requiredInput)) {
        st <- st[grepl(
            paste(requiredInput, collapse = "|"), st$requiredInput), ]
    }
    if (!description) {
        st <- st[, -which(colnames(st)=="description")]
    }
    return(st)
}


#' Multiple Signatures Computation
#'
#' @description This function computes all the signatures for a specific
#' 'inputType'. Further, it is possible to select specific signatures
#' setting the 'tumor', the 'tissue' and/or the 'topic'.
#'
#' @param dataset Expression values. A data frame or a matrix where rows
#' correspond to genes and columns correspond to samples.
#' Alternatively, an object of type \linkS4class{SummarizedExperiment},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}.
#' @param nametype character string saying the type of gene name ID (row names
#' in dataset). Either one of "SYMBOL", "ENTREZID" or "ENSEMBL".
#' @param inputType character vector saying the type of data you are using.
#' When working with bulk data this should be either one of "microarray" or
#' "rnaseq". When working with single-cell data and spatial transcriptomics data
#' this could be "sc" to compute only signatures developed by single-cell data
#' or c("rnaseq", "sc") to compute all the signatures.
#' @param whichAssay integer scalar or string indicating which assay of
#' dataset to use.
#' @param whichSign character vector saying the signatures to compute.
#' @param tumor character vector saying the tumor types. Signatures from that
#' tumors will be computed (this can also be "pan-cancer").
#' @param tissue character vector saying the tumor tissues. Signatures from
#' that tissues will be computed (this can also be "pan-tissue").
#' @param topic character vector saying signatures topics. Signatures having
#' that topics will be computed.
#' @param ... other arguments passed on to the signature functions.
#'
#' @return A SummarizedExperiment object in which the signatures' scores
#' are added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' data(ovse)
#' multipleSign(dataset = ovse)
#' multipleSign(dataset = ovse, tissue = "ovary")
#'
#' @export
multipleSign <- function(
        dataset, nametype = "SYMBOL", inputType = "rnaseq",
        whichAssay = "norm_expr", whichSign = NULL, tumor = NULL,
        tissue = NULL, topic = NULL, ...) {
    st <- availableSignatures(
        tumor = tumor,
        tissue = tissue,
        topic = topic,
        requiredInput = inputType)

    inputType <- inputType[!inputType=="sc"]
    argg <- c(as.list(environment()), list(...))
    argg <- argg[!names(argg)=="st"]

    if (!is.null(whichSign)) {
        st <- st[grepl(paste(whichSign, collapse = "|"), st$signature), ]
    }

    for (i in seq_len(nrow(st))) {
        author <- st[i, "author"]
        argg_i <- c(argg, list(author = author))

        signName <- st[i, "functionName"]
        signFunc <- eval(parse(text = signName))
        signFuncArg <- formals(signFunc)

        argg_i[!(names(argg_i) %in% names(signFuncArg))] <- NULL
        signFuncArg[names(signFuncArg) %in% names(argg_i)] <- NULL
        matchedArg <- c(argg_i, signFuncArg)
        matchedArg[names(matchedArg) == "..."] <- NULL

        dataset <- tryCatch(
            do.call(signFunc, matchedArg),
            error = function(e){
                message(signName, author, " raised an error:\n",
                        e$message,
                        "\nTherefore it was omitted")
                return(dataset)
            })
        argg$dataset <- dataset
    }
    return(dataset)
}


#' Get Signature Gene List
#'
#' @description This function returns the list of genes of a signature.
#'
#' @param whichSign name of the signature. The names are those in column
#' 'signature' from the table which is obtained by
#' \code{\link[signifinder]{availableSignatures}}.
#'
#' @return A dataframe object with "SYMBOL" in the first column. Some signatures
#' have also additional colums: "coeff" for coefficients that weigh the gene
#' contributions; "class" for a classification that divides the signature in
#' two or more groups. Few signatures have other specific columns.
#'
#' @examples
#' getSignGenes("EMT_Miow")
#'
#' @export
getSignGenes <- function(whichSign) {
    if (!(whichSign %in% signatureTable$signature)) {
        stop("This signature is not present in signifinder")}

    if(whichSign == "Combined_Thompson"){
        message(
        "Combined_Thompson is a combination of EMT_Thompson and
        Tinflam_Thompson signatures.")
    } else {
        signGene <- eval(parse(text = whichSign))
        rownames(signGene) <- seq_len(nrow(signGene))
        return(signGene)
    }
}

#' Pancancer Cellular States Barkley Scoring Function
#'
#' @importFrom ggplot2 cut_number
#' @export
.barkleyFun <- function(dataset, signList, modules) {
  
  data.avg <- sort(rowMeans(x = dataset, na.rm = TRUE))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = 25, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  mod <- signList[modules]
  
  rand = lapply(names(mod), function(m){
    lapply(1:1000, function(i){
      used = vector()
      unused = binned
      for (g in mod[[m]]){
        pool = data.cut[g]
        if (!(is.na(pool))) {
          new = sample(unused[[pool]], 1)
          used = c(used, new)
          unused[[pool]] = setdiff(unused[[pool]], new)
        }
      }
      used})
  })
  names(rand) = names(mod)
  
  scores = t(sapply(names(mod), function(m){
    ra = sapply(rand[[m]], function(i){
      colMeans(dataset[i, ], na.rm = TRUE)
    })
    re = colMeans(dataset[rownames(dataset) %in% mod[[m]], ], na.rm = TRUE)
    p = rowMeans(ra >= re)
    p = -log10(p)
  }))
  scores[is.infinite(scores)] = 4
  scores = scores/4
  
  return(scores)
}
