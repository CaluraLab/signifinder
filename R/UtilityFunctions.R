SignatureNames <- c(
    "EMT_Miow_Epithelial",
    "EMT_Miow_Mesenchymal",
    "EMT_Mak",
    "EMT_Cheng",
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
    "TLS_Cabrita",
    "StemCellCD49f_Smith",
    "CIN_Carter_25",
    "CIN_Carter_70",
    "CellCycle_Lundberg",
    "CellCycle_Davoli",
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
    "IPSOV"
)

#' @importFrom viridis magma viridis
mycol <- c("#FCFDD4", rev(viridis::magma(10)))
mycol1 <- rev(viridis::viridis(10))

my_colors <- RColorBrewer::brewer.pal(5, "Spectral")

my_colors <- colorRampPalette(my_colors)(100)

.GetGenes <- function(name) {
    if (name %in% c("EMT_Miow_Epithelial", "EMT_Miow_Mesenchymal")) {
        sname <- paste0(substring(name, 10), "-like")
        g <- EMT_Miow$SYMBOL[EMT_Miow$class == sname]
    } else if (name %in% c("ECM_Chakravarthy_up", "ECM_Chakravarthy_down")) {
        sname <- substring(name, 18)
        g <- ECM_Chakravarthy$SYMBOL[grepl(sname, ECM_Chakravarthy$class)]
    } else if (name %in% c(
        "ConsensusOV_Chen_IMR",
        "ConsensusOV_Chen_DIF",
        "ConsensusOV_Chen_PRO",
        "ConsensusOV_Chen_MES"
    )) {
        sname <- substring(name, 18)
        g <- ConsensusOV_Chen$SYMBOL[ConsensusOV_Chen$class == sname]
    } else if (name %in% c(
        "IPS_Charoentong_MHC",
        "IPS_Charoentong_CP",
        "IPS_Charoentong_EC",
        "IPS_Charoentong_SC"
    )) {
        g <-IPS_Charoentong$SYMBOL[IPS_Charoentong$class == substring(name, 17)]
    } else if (name %in% c("CIN_Carter_25", "CIN_Carter_70")) {
        g <- if (name == "CIN_Carter_25") {
            CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"]
        } else {
            CIN_Carter$SYMBOL
        }
    } else if (name %in% c(
        "ISC_MerlosSuarez_ISCEphB2",
        "ISC_MerlosSuarez_LateTA",
        "ISC_MerlosSuarez_ISCLgr5",
        "ISC_MerlosSuarez_Prolif"
    )) {
        sname <- substring(name, 18)
        g <- ISC_MerlosSuarez$SYMBOL[ISC_MerlosSuarez$class == sname]
    } else if (name == "Tinflam_Ayers") {
        g <- Tinflam_Ayers$SYMBOL[Tinflam_Ayers$class == "TInflam"]
    } else if (name %in% c("Autophagy_ChenM_OS", "Autophagy_ChenM_DFS")) {
        sname <- substring(name, 17)
        g <- Autophagy_ChenM$SYMBOL[Autophagy_ChenM$class == sname]
    } else if (name %in% c(
        "EMT_Mak",
        "EMT_Cheng",
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
        "Matrisome_Yuzhalin",
        "Chemokines_Messina",
        "ImmunoScore_Hao",
        "ImmunoScore_Roh",
        "ImmuneCyt_Rooney",
        "ImmuneCyt_Davoli",
        "IFN_Ayers",
        "HRDS_Lu",
        "ExpandedImmune_Ayers",
        "CellCycle_Davoli",
        "PassON_Du",
        "VEGF_Hu",
        "DNArep_Kang",
        "ASC_Smith",
        "IPS_Charoentong",
        "StemCellCD49f_Smith",
        "Glycolysis_Zhang",
        "Glycolysis_Xu",
        "MitoticIndex_Yang",
        "Autophagy_Xu",
        "Autophagy_Wang",
        "Autophagy_ChenH",
        "CellCycle_Lundberg",
        "IPRES_Hugo",
        "CIS_Robertson",
        "TLS_Cabrita",
        "IPSOV_Shen"
    )) {
        g <- unique(eval(parse(text = name))[, "SYMBOL"])
    }
    res <- cbind(g, rep(name, length(g)))
    colnames(res) <- c("Gene", "Signature")
    return(res)
}

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
    if (!all(sName %in% SignatureNames)) {
        stop(paste(
            "signatures must be among:",
            paste(SignatureNames, collapse = ", ")
        ))
    }
    if (!all(sName %in% colnames(colData(data)))) {
        stop("signature names must be in data")
    }
}

.matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
.getMatrix <- function(userdata) {
    if (!is.matrix(userdata)) {
        if (class(userdata) %in% c(
            "SpatialExperiment",
            "SummarizedExperiment",
            "SingleCellExperiment")) {
            userdata <- as.matrix(assay(userdata))
        } else if (is.data.frame(userdata)) {
            userdata <- as.matrix(userdata)
        } else { stop("This dataset type is not supported") }
    }
    return(userdata)
}

#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<-
.returnAsInput <- function(userdata, result, SignName, datasetm) {
    if (!is.matrix(userdata) & !is.data.frame(userdata)) {
        if (class(userdata) %in% c(
            "SpatialExperiment",
            "SummarizedExperiment",
            "SingleCellExperiment")) {
            names <- c(colnames(userdata@colData), SignName)
            if (is.vector(result)) {
                userdata@colData <- cbind(userdata@colData, name = result)
                colnames(userdata@colData) <- names
            } else {
                userdata@colData <- cbind(userdata@colData, t(result))
            }
        }
        return(userdata)
    } else if (is.matrix(userdata) | is.data.frame(userdata)) {
        if (is.vector(result)) {
            result <- SummarizedExperiment(
                assays = list(norm_expr = datasetm),
                colData = data.frame(name = result)
            )
            colnames(colData(result)) <- SignName
        } else {
            result <- SummarizedExperiment(
                assays = list(norm_expr = datasetm),
                colData = t(result)
            )
        }
        return(result)
    }
}

.ipsmap <- function(x) {
    if (is.na(x)) {
        NA
    } else {
        if (x <= 0) {
            0
        } else if (x >= 3) {
            10
        } else {
            round(x*10/3, digits = 0) }}
}

#' @importFrom parallel mclapply
.GSVAPvalues <- function(expr, gset.idx.list, gsvaResult, nperm, args) {
    datasetGenes <- rownames(expr)
    filteredGeneSets <- lapply(gset.idx.list, y = datasetGenes, intersect)
    permutedResults <- mclapply(seq_len(nperm), function(x) {
        message("Performing permutation number", x, "\n")
        permlist <- lapply(seq_len(length(gset.idx.list)), function(i) {
            sample(
                datasetGenes,
                size = lengths(filteredGeneSets)[i],
                replace = FALSE) })
        args$gset.idx.list <- permlist
        gsva_matrix <- do.call(gsva, args)
        data.frame(t(gsva_matrix))
    }, mc.cores = 1)
    permutedResByGeneSet <- split.default(x = Reduce(
        cbind, permutedResults), seq_len(length(gset.idx.list)))
    permutedResByGeneSet <- lapply(permutedResByGeneSet, function(x) {
        data.frame(t(x)) })
    finalRes <- do.call(
        rbind, lapply(seq_len(length(gset.idx.list)), function(i) {
            gspvalues <- sapply(seq_len(ncol(expr)), function(j) {
                (min(c(
                    sum(permutedResByGeneSet[[i]][, j] <= gsvaResult[i, j]),
                    sum(permutedResByGeneSet[[i]][, j] >= gsvaResult[i, j])
                )) + 1) / (nperm + 1) })
            gspvalues })
    )
    colnames(finalRes) <- colnames(expr)
    rownames(finalRes) <- paste(names(gset.idx.list), "pval", sep = "_")
    return(finalRes)
}

.consistencyCheck <- function(nametype, functionName, author = NULL) {
    if (!(nametype %in% c("SYMBOL", "ENTREZID", "ENSEMBL"))) {
        stop("The name of genes must be either SYMBOL, ENTREZID or ENSEMBL")}
    if (!is.null(author)) {
        if (!(author %in% signatureTable$author[
            signatureTable$functionName==functionName])) {
            stop("This author is not present for this signature")}}
}

.coeffScore <- function(
        sdata, datasetm, namesignature, detail = NULL, author = "") {
    .percentageOfGenesUsed(
        namesignature, datasetm, sdata$SYMBOL, detail, author)
    sdata <- sdata[sdata$SYMBOL %in% row.names(datasetm), ]
    columnNA <- .managena(datasetm = datasetm, genes = sdata$SYMBOL)
    score <- colSums(datasetm[sdata$SYMBOL, ] * sdata$coeff, na.rm = TRUE)
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
    columnNA <- (length(genes) - colSums(!is.na(datasetm))) / length(genes)
    if (sum(columnNA > 0.9) > 0) {
        warning("Some samples have more than 90% NA values")}
    return(columnNA)
}

#' @importFrom SummarizedExperiment assays SummarizedExperiment assays<-
#' @importFrom BiocGenerics width
#' @importFrom IRanges reduce
#' @importFrom DGEobj.utils convertCounts
.dataTransformation <- function(dataset, data, trans, hgReference, nametype) {
    if (class(dataset)[1] %in% c(
        "SpatialExperiment",
        "SummarizedExperiment",
        "SingleCellExperiment"
    )) { if (trans %in% names(assays(dataset))) { return(dataset) } }

    if (hgReference == "hg19") {
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (hgReference == "hg38") {
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else { stop("Human reference genome must be either hg19 or hg38") }

    exons.db <- ensembldb::exonsBy(txdb, by = "gene")

    g <- rownames(data)
    egs <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, keys = g, column = "ENTREZID",
        keytype = nametype, multiVals = "first")
    use <- !is.na(egs)

    exons_g <- NULL
    exons_g[use] <- lapply(egs[use], function(eg) { exons.db[[eg]] })
    names(exons_g) <- names(egs)
    use <- !sapply(exons_g, is.null)

    glen <- NULL
    glen[use] <- sapply(names(egs[use]), function(eg) {
        sum(width(reduce(exons_g[[eg]])))
    })

    usec <- colSums(data[use,])>0

    tdata <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    tdata[use,usec] <- convertCounts(
        countsMatrix = data[use,usec],
        unit = trans, geneLength = glen[use])

    if (class(dataset)[1] %in% c(
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
        name, author, " is using ", round(g_per), "% of genes")
        if (g_per < 30) { warning(
            name, author, " is computed with less than 30% of its genes")}
    } else { message(
        name, author, " is using ", round(g_per), "% of ", detail, " genes")
        if (g_per < 30) { warning(
            name, author, detail,
            " is computed with less than 30% of its genes")}}
    # return(g_per==0)
}

.geneIDtrans <- function(nametype, genes) {
    if (nametype != "SYMBOL") { AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, keys = genes, column = nametype,
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
        rmatrix <- matrix(rnorm(n, 100, 50), ncol = 5)
    }
    rownames(rmatrix) <- g
    return(rmatrix)
}


#' Show Available Signatures
#'
#' It shows a table containing all the information of the signatures collected
#' in the package.
#'
#' @param tumor character vector saying the type of tumors for which signatures
#' are developed. Used to filter the signatures in the table.
#' @param tissue character vector saying the type of tissues for which signatures
#' are developed. Used to filter the signatures in the table.
#' @param topic character vector saying the signature topics. Used to filter
#' the signatures in the table.
#' @param requiredInput character string saying the type of data required in
#' input by the signature. Either one of "microarray" or "rnaseq". Used to
#' filter the signatures in the table.
#' @param description logical. If TRUE it shows the signature's description.
#'
#' @return A data frame with 46 rows and 10 variables:
#' \describe{
#'   \item{signature}{name of the signature}
#'   \item{scoreLabel}{label of the signature when added inside colData section}
#'   \item{functionName}{name of the function to use to compute the signature}
#'   \item{topic}{main cancer topic of the signature}
#'   \item{tumor}{tumor type for which the signature was developed}
#'   \item{tissue}{tumor tissue for which the signature was developed}
#'   \item{requiredInput}{type of data with which the signature was developed}
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
        st <- st[, -10]
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
#' @param inputType character string saying the type of data you are using.
#' Either one of "microarray" or "rnaseq".
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
#' multipleSign(dataset = ovse)
#' multipleSign(dataset = ovse, tissue = "ovary")
#'
#' @export
multipleSign <- function(
        dataset, nametype = "SYMBOL", inputType = "rnaseq",
        whichSign = NULL, tumor = NULL, tissue = NULL, topic = NULL, ...) {
    argg <- c(as.list(environment()), list(...))
    st <- availableSignatures(
        tumor = tumor,
        tissue = tissue,
        topic = topic,
        requiredInput = inputType
    )

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
