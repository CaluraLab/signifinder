
#' Epithelial-Mesenchymal Transition Signature
#'
#' @description This signature is computed accordingly to the reference paper,
#' to have more details explore the function
#'  \code{\link[signifinder]{availableSignatures}}.
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
#' @param tumorTissue type of tissue for which the signature is developed.
#' @param author first author of the specific signature publication.
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the Epithelial and Mesenchymal
#' scores are added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
EMTSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray",
    tumorTissue = "ovary", author = "Miow", pvalues = FALSE, nperm = 100, ...) {

    firstCheck(nametype, tumorTissue, "EMTSign", author = author)

    datasetm <- getMatrix(dataset)

    if(tumorTissue == "ovary" & author == "Miow"){

       if(nametype!="SYMBOL"){
        EMTMiowdata$Gene_Symbol <- mapIds(
            org.Hs.eg.db, keys = EMTMiowdata$Gene_Symbol, column = nametype,
            keytype = "SYMBOL", multiVals = "first")}

        Signature_EL <- EMTMiowdata[grep('Epithelial', EMTMiowdata$Category),]
        Signature_ML <- EMTMiowdata[grep('Mesenchymal', EMTMiowdata$Category),]

        percentageOfGenesUsed("EMTSign", datasetm, Signature_EL$Gene_Symbol,
                                detail = "epithelial-like")
        percentageOfGenesUsed("EMTSign", datasetm, Signature_ML$Gene_Symbol,
                                detail = "mesenchymal-like")

        gene_sets <- list(Epithelial = Signature_EL$Gene_Symbol,
                        Mesenchymal = Signature_ML$Gene_Symbol)

        dots <- list(...)
        kcdftype <- ifelse(inputType == "microarray", "Gaussian", "Poisson")
        args <- matchArguments(dots, list(expr = datasetm,
                                        gset.idx.list = gene_sets,
                                        method = "ssgsea", kcdf = kcdftype,
                                        min.sz = 5, ssgsea.norm = FALSE,
                                        verbose = FALSE))
        gsva_matrix <- suppressWarnings(do.call(gsva, args))

        if(pvalues){
            gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets,
                                    gsvaResult = gsva_matrix,
                                    nperm = nperm, args = args)
            gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

        return(returnAsInput(userdata = dataset, result = gsva_matrix,
                            SignName = "", datasetm))

    } else {
        if(tumorTissue == "pan-tissue" & author == "Mak") {
            if(nametype!="SYMBOL"){
                EMTMakdata$Gene_Symbol <- mapIds(
                    org.Hs.eg.db, keys = EMTMakdata$Gene_Symbol, column = nametype,
                    keytype = "SYMBOL", multiVals = "first")}

            Sign_E <- EMTMakdata[grep('E', EMTMakdata$Category),]
            Sign_M <- EMTMakdata[grep('M', EMTMakdata$Category),]

            percentageOfGenesUsed("EMTSign", datasetm, Sign_E$Gene_Symbol,
                                  detail = "epithelial-like")
            percentageOfGenesUsed("EMTSign", datasetm, Sign_M$Gene_Symbol,
                                  detail = "mesenchymal-like")

            datasetm_n <- datasetm[intersect(
                row.names(datasetm), EMTMakdata$Gene_Symbol), ]
            datasetm_n <- if(inputType == "rnaseq") {log2(datasetm_n)
                } else {datasetm_n}
            columnNA <- managena(datasetm_n, genes = EMTMakdata$Gene_Symbol)
            score <- colMeans(
                datasetm_n[intersect(Sign_M$Gene_Symbol, row.names(datasetm_n)),]
                ) - colMeans(
                datasetm_n[intersect(Sign_E$Gene_Symbol, row.names(datasetm_n)), ])
            score[columnNA > 0.9] <- NA

        } else if(tumorTissue == "breast" & author == "Cheng") {
            if(nametype!="SYMBOL") {
                EMTChengdata <- mapIds(
                    org.Hs.eg.db, keys = EMTChengdata, column = nametype,
                    keytype = "SYMBOL", multiVals = "first")}

            percentageOfGenesUsed("EMTSign", datasetm, EMTChengdata)

            datasetm_n <- datasetm[intersect(row.names(datasetm), EMTChengdata), ]
            datasetm_n <- if(inputType == "rnaseq") {log2(datasetm_n)
                } else {datasetm_n}
            columnNA <- managena(datasetm = datasetm_n, genes = EMTChengdata)
            score <- prcomp(t(datasetm_n))$x[,1]
            score[columnNA > 0.9] <- NA
        }
        return(returnAsInput(userdata = dataset, result = score,
                             SignName = paste0("EMT", author), datasetm))}
}


#' Pyroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the Pyroptosis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
pyroptosisSign <- function(dataset, nametype = "SYMBOL", inputType = "rnaseq",
                tumorTissue = "ovary", author = "Ye", hgReference = "hg19"){

    firstCheck(nametype, tumorTissue, "pyroptosisSign", author)

    Pyroptosisdata <- get(paste0("Pyroptosis", author, "data"))

    datasetm <- getMatrix(dataset)

    if(tumorTissue == "ovary" & author == "Ye"){
        # datasetm_n <- scale(datasetm)
        datasetm_n <- dataTransformation(datasetm, "FPKM", hgReference)
    } else if (tumorTissue == "stomach" & author == "Shao"){
        if(inputType == "rnaseq"){
            datasetm_n <- dataTransformation(datasetm, "FPKM", hgReference)}
    } else if (tumorTissue == "lung" & author == "Lin"){
        datasetm_n <- dataTransformation(datasetm, "TPM", hgReference)
    } else {datasetm_n <- datasetm}

    score <- coefficientsScore(Pyroptosisdata, datasetm_n, nametype,
                                "pyroptosisSign")

    return(returnAsInput(
        userdata = dataset, result = score,
        SignName = paste0("Pyroptosis", author), datasetm))
}


#' Ferroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Ferroptosis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ferroptosisSign <- function(dataset, nametype = "SYMBOL",
                            tumorTissue = "ovary", author = "Ye"){

    firstCheck(nametype, tumorTissue, "ferroptosisSign", author)
    Ferroptosisdata <- get(paste0("Ferroptosis", author, "data"))

    datasetm <- getMatrix(dataset)
    score <- coefficientsScore(Ferroptosisdata, datasetm,
                                nametype, "ferroptosisSign")

    if(tumorTissue == "liver" & author == "Liang" ){score <- exp(score)}

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = paste0("Ferroptosis", author), datasetm))
}


#' Lipid Metabolism Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Lipid scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
lipidMetabolismSign <- function(dataset, nametype = "SYMBOL",
                                tumorTissue = "ovary") {

    firstCheck(nametype, tumorTissue, "lipidMetabolismSign")

    datasetm <- getMatrix(dataset)
    score <- coefficientsScore(LipidMetabolismdata, datasetm,
                                nametype, "lipidMetabolismSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "LipidMetabolism", datasetm))
}


#' Hypoxia Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Hypoxia scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom matrixStats colMedians
#' @importFrom stats setNames
#'
#' @import org.Hs.eg.db
#'
#' @export
hypoxiaSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray",
                        tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "hypoxiaSign")

    if(nametype=="SYMBOL") {genetouse <- Hypoxiadata$Gene_Symbol
    } else if(nametype=="ENSEMBL") {genetouse <- Hypoxiadata$Gene_Ensembl
    } else {genetouse <- mapIds(org.Hs.eg.db, keys = Hypoxiadata$Gene_Symbol,
        column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)
    datasetm_n <- if(inputType == "rnaseq") {log2(datasetm)} else {datasetm}
    score <- statScore(genetouse, datasetm = datasetm_n, nametype = "SYMBOL",
        typeofstat = "median", namesignature = "hypoxiaSign")

    return(returnAsInput(userdata = dataset, result = as.vector(scale(score)),
        SignName = "Hypoxia", datasetm))
}


#' Platinum Resistance Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the Platinum Resistance scores
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
platinumResistanceSign <- function(dataset, nametype = "SYMBOL",
        tumorTissue = "ovary", pvalues = FALSE, nperm = 100, ...){

    firstCheck(nametype, tumorTissue, "platinumResistanceSign")

    if(nametype!="SYMBOL"){
        PlatinumResistancedata$Gene_Symbol <- mapIds(
            org.Hs.eg.db, keys = PlatinumResistancedata$Gene_Symbol,
            column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    Signature_up <- PlatinumResistancedata[
        grep('PlatinumResistanceUp', PlatinumResistancedata$Category),]
    Signature_down <- PlatinumResistancedata[
        -grep('PlatinumResistanceUp', PlatinumResistancedata$Category),]

    percentageOfGenesUsed("platinumResistanceSign", datasetm,
                          Signature_up$Gene_Symbol, detail = "up")
    percentageOfGenesUsed("platinumResistanceSign", datasetm,
                          Signature_down$Gene_Symbol, detail = "down")

    gene_sets <- list(PlatinumResistanceUp = Signature_up$Gene_Symbol,
                      PlatinumResistanceDown = Signature_down$Gene_Symbol)

    dots <- list(...)
    args <- matchArguments(dots, list(
        expr = datasetm, gset.idx.list = gene_sets,
        method = "gsva", kcdf = "Gaussian", min.sz = 5,
        ssgsea.norm = FALSE, verbose = FALSE))
    gsva_count <- suppressWarnings(do.call(gsva, args))

    if(pvalues){
        gsva_pval <- GSVAPvalues(
            expr = datasetm, gset.idx.list = gene_sets,
            gsvaResult = gsva_matrix, nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(
        userdata = dataset, result = gsva_count, SignName = "", datasetm))
}


#' Immunogenic Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Immunogenic scores will
#' be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
immunoScoreSign <- function(dataset, nametype = "SYMBOL",
                            tumorTissue = "ovary", author = "Hao"){

    firstCheck(nametype, tumorTissue, "immunoScoreSign", author)

    datasetm <- getMatrix(dataset)

    if(tumorTissue == "ovary"){
        if(nametype!="SYMBOL"){
            ImmunoScoreHaodata$genes <- mapIds(
                org.Hs.eg.db, keys = ImmunoScoreHaodata$genes,
                column = nametype, keytype = "SYMBOL", multiVals = "first")}

        g <- intersect(row.names(datasetm), ImmunoScoreHaodata$genes)

        percentageOfGenesUsed("immunoScoreSign", datasetm,
                              ImmunoScoreHaodata$genes)

        datasetm_n <- datasetm[g,]
        ImmunoScoreHaodata <- ImmunoScoreHaodata[ImmunoScoreHaodata$genes%in%g,]
        columnNA <- managena(datasetm_n, g)
        SE <- (ImmunoScoreHaodata$HR - ImmunoScoreHaodata$`95CI_L`)/1.96
        k <- (1 - ImmunoScoreHaodata$HR)/SE
        score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(p){
                sum(k*datasetm_n[,p], na.rm = TRUE)}))
        score[columnNA > 0.9] <- NA
    } else if(tumorTissue=="pan-tissue"){
        score <- statScore(
            ImmunoScoreRohdata, datasetm = datasetm, nametype = nametype,
            typeofstat = "meang", namesignature = "immunoScoreSign")}

    return(returnAsInput(userdata = dataset, result = score,
        SignName = paste0("ImmunoScore", author), datasetm))
}


#' ConsensusOV Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param method the subtyping method to use. Default is "consensusOV".
#' @param ... optional parameters to be passed to
#' \code{\link[consensusOV]{get.subtypes}}.
#'
#' @return A SummarizedExperiment object in which the COnsensusOV scores
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom consensusOV get.subtypes
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
consensusOVSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary",
                            method = "consensusOV", ...){

    firstCheck(nametype, tumorTissue, "consensusOVSign")

    datasetm <- getMatrix(dataset)

    if(nametype!="ENTREZID"){
        genename <- mapIds(org.Hs.eg.db, keys = row.names(datasetm),
                            column = "ENTREZID", keytype = nametype,
                            multiVals = "first")
        datasetm_n <- datasetm[!is.na(genename),]
        genename <- genename[!is.na(genename)]
        datasetm_n <- datasetm_n[!duplicated(genename),]
        genename <- genename[!duplicated(genename)]
    } else {
        genename <- row.names(datasetm)
        datasetm_n <- datasetm}

    consensus_subtypes <- get.subtypes(expression.dataset = datasetm_n,
                            entrez.ids = genename, method = method, ...)

    return(returnAsInput(userdata = dataset,
                        result = t(consensus_subtypes$rf.probs),
                        SignName = "", datasetm))
}


#' ImmunoPhenoScore Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the
#' IPS, MHC, CP, EC and SC scores will be added in the
#' \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices pdf dev.off
#' @importFrom stats sd
#'
#' @export
IPSSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "IPSSign")

    if(nametype!="SYMBOL"){
        IPSdata[,c(1,2)] <- data.frame(lapply(IPSdata[,c(1,2)], function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype,
                                    keytype="SYMBOL", multiVals="first"))))}

    datasetm <- getMatrix(dataset)
    sample_names <- colnames(datasetm)

    percentageOfGenesUsed("IPSSign", datasetm, IPSdata$GENE)

    ind <- which(is.na(match(IPSdata$GENE, row.names(datasetm))))
    MISSING_GENES <- IPSdata$GENE[ind]
    if (length(MISSING_GENES)>0) {cat("Differently named or missing genes: ",
                                    MISSING_GENES, "\n")}

    IPS <- NULL; MHC <- NULL; CP <- NULL; EC <- NULL; SC <- NULL; AZ <- NULL
    for (i in 1:length(sample_names)) {
        GE <- datasetm[,i]
        Z1 <- (datasetm[match(
            IPSdata$GENE, row.names(datasetm)),i]-mean(GE))/sd(GE)
        WEIGHT <- NULL; MIG <- NULL; k <- 1
        for (gen in unique(IPSdata$NAME)) {
            MIG[k] <- mean(Z1[which(IPSdata$NAME==gen)], na.rm=TRUE)
            WEIGHT[k] <- mean(IPSdata$WEIGHT[which(IPSdata$NAME==gen)])
            k<-k+1}
        WG <- MIG*WEIGHT
        MHC[i]<-mean(WG[1:10], na.rm = TRUE)
        CP[i]<-mean(WG[11:20], na.rm = TRUE)
        EC[i]<-mean(WG[21:24], na.rm = TRUE)
        SC[i]<-mean(WG[25:26], na.rm = TRUE)
        AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
        IPS[i]<-ipsmap(AZ[i])}

    ipsres <- data.frame(IPS, MHC, CP, EC, SC)
    row.names(ipsres) <- sample_names
    return(returnAsInput(userdata = dataset, result = t(ipsres),
                        SignName = "", datasetm))
}


#' Core Matrisome Gene signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom matrixStats colMedians
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
matrisomeSign <- function(dataset, nametype = "SYMBOL",
                        tumorTissue = "pan-tissue") {

    firstCheck(nametype, tumorTissue, "matrisomeSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(Matrisomedata, datasetm, nametype,
                        "median", "matrisomeSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "Matrisome", datasetm))
}


#' Mitotic Index
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
mitoticIndexSign <- function(dataset, nametype = "SYMBOL",
                            tumorTissue = "pan-tissue") {

    firstCheck(nametype, tumorTissue, "mitoticIndexSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(MitoticIndexdata, datasetm, nametype,
                        "mean", "mitoticIndexSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "MitoticIndex", datasetm))
}


#' Immune Cytolytic Activity Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
ImmuneCytSign <- function(dataset, nametype = "SYMBOL",
                          tumorTissue = "pan-tissue", author = "Rooney"){

    firstCheck(nametype, tumorTissue, "ImmuneCytSign", author)

    datasetm <- getMatrix(dataset)
    if(tumorTissue == "pan-tissue" & author == "Rooney"){
         score <- statScore(
             ImmuneCytRooneydata, datasetm, nametype, "meang", "ImmuneCytSign")
    } else if(tumorTissue == "pan-tissue" & author == "Davoli") {
        score <- statScore(
            ImmuneCytDavolidata, datasetm, nametype, "mean", "ImmuneCytSign")
    }
    return(returnAsInput(userdata = dataset, result = score,
                         SignName = paste0("ImmuneCyt", author), datasetm))
}


#' IFN-γ Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
IFNSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "IFNSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(IFNdata, datasetm, nametype, "mean", "IFNSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "IFN", datasetm))
}


#' ExpandedImmune Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
expandedImmuneSign <- function(dataset, nametype = "SYMBOL",
                                tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "expandedImmuneSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(ExpandedImmunedata, datasetm, nametype,
                        "mean", "expandedImmuneSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "ExpandedImmune", datasetm))
}


#' Tertiary Lymphoid Structures (TLS) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
TLSSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "skin"){

    firstCheck(nametype, tumorTissue, "TLSSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(TLSdata, datasetm, nametype, "meang", "TLSSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "TLS", datasetm))
}


#' CD49fHi Basal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
CD49BSCSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "prostate"){

    firstCheck(nametype, tumorTissue, "CD49BSCSign")

    datasetm <- getMatrix(dataset)
    score <- coefficientsScore(CD49BSCdata, datasetm, nametype, "CD49BSCSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "CD49BSC", datasetm))
}


#' Chromosomal instability Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param typeofCIN A pre-defined number 25 or 70. It represents the number of
#' genes on which the score can be calculated. The 25 genes represent a subgroup
#' derived from the 70.
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
CINSign <- function(dataset, nametype = "SYMBOL",
                    tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "CINSign")

    datasetm <- getMatrix(dataset)

    score25 <- statScore(CINdata$SYMBOL[CINdata$class == "CIN25"],
                        scale(datasetm, center = TRUE, scale = FALSE),
                        nametype, "sum", "CINSign")
    score70 <- statScore(CINdata$SYMBOL,
                        scale(datasetm, center = TRUE, scale = FALSE),
                        nametype, "sum", "CINSign")
    score <- data.frame(CIN25=score25, CIN70=score70)

    return(returnAsInput(userdata = dataset, result = t(score),
                        SignName = "", datasetm))
}


#' Cell-cycle Signature classifier
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
CCSSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue",
                    author = "Lundberg"){

    firstCheck(nametype, tumorTissue, "CCSSign", author)

    datasetm <- getMatrix(dataset)

    if(author == "Lundberg"){
        score <- statScore(CCSLundbergdata, datasetm, nametype,
                            "sum", "CCSSign")
    } else if(author == "Davoli"){
        score <- statScore(CCSDavolidata, datasetm, nametype,
                            "mean", "CCSSign")}

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = paste0("CCS", author), datasetm))
}


#' Chemokine Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
chemokineSign <- function(dataset, nametype = "SYMBOL",
                          tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "chemokineSign")

    if(nametype!="SYMBOL") {Chemokinedata <- mapIds(
        org.Hs.eg.db, keys = Chemokinedata, column = nametype,
        keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ChemokineSign", datasetm, Chemokinedata)

    datasetm_n <- datasetm[intersect(row.names(datasetm), Chemokinedata), ]
    columnNA <- managena(datasetm_n, Chemokinedata)
    score <- prcomp(t(datasetm), center = TRUE, scale = TRUE)$x[, 1]
    score[columnNA > 0.9] <- NA

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "Chemokine", datasetm))
}

#' Adult Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ASCSign <- function(dataset, nametype= "SYMBOL",
                    tumorTissue = "epithelial-derived neuroendocrine cancer"){

    firstCheck(nametype, tumorTissue, "ASCSign")

    if(nametype!="SYMBOL"){
        ASCdata <- mapIds(org.Hs.eg.db, keys = ASCdata, column = nametype,
                          keytype = "SYMBOL", multiVals = "first")}
    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ASCSign", datasetm, ASCdata)

    datasetm_n <- log2(datasetm[row.names(datasetm) %in% ASCdata, ] + 1)
    columnNA <- managena(datasetm_n, ASCdata)
    score <- colSums((datasetm_n - rowMeans(datasetm_n))/
                    sapply(as.data.frame(t(datasetm_n)), sd, na.rm = TRUE))
    score[columnNA > 0.9] <- NA

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = "ASC", datasetm))
}


#' passON Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the passON score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
PassONSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "skin", ...){

    firstCheck(nametype, tumorTissue, "PassONSign")

    if(nametype!="SYMBOL"){
        PASS.ONdata <- lapply(PASS.ONdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype,
                                    keytype="SYMBOL", multiVals="first")))}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("PassONSign", datasetm, unlist(PASS.ONdata))

    dots <- list(...)
    args <- matchArguments(dots, list(
        expr = datasetm, gset.idx.list = PASS.ONdata, method = "ssgsea",
        kcdf = "Poisson", min.sz = 5, ssgsea.norm = TRUE, verbose = FALSE))

    gsva_matrix <- suppressWarnings(do.call(gsva, args))

    gsva_mean <- sapply(colnames(gsva_matrix), function(x) {
        weighted.mean(gsva_matrix[,x], c(23, 26, 81, 18))
    })

    return(returnAsInput(userdata = dataset, result = gsva_mean,
                            SignName = "PASS.ON", datasetm))
}

#' IPRES Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the passON score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
IPRESSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "skin",
                      pvalues = FALSE, nperm = 100, ...) {

    firstCheck(nametype, tumorTissue, "IPRESSign")

    if(nametype!="SYMBOL"){
        IPRESdata <- lapply(IPRESdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype,
                                    keytype="SYMBOL", multiVals="first")))}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("IPRESSign", datasetm, unlist(ipresdata))

    dots <- list(...)
    args <- matchArguments(dots, list(
        expr = datasetm, gset.idx.list = IPRESdata, method = "ssgsea",
        kcdf = "Poisson", min.sz = 5, ssgsea.norm = TRUE, verbose = FALSE))

    gsva_matrix <- suppressWarnings(do.call(gsva, args))

    if(pvalues){
        gsva_pval <- GSVAPvalues(
            expr = datasetm, gset.idx.list = gene_sets,
            gsvaResult = gsva_matrix, nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    z_score <- (gsva_matrix - rowMeans(gsva_matrix))/sapply(
        as.data.frame(t(gsva_matrix)), sd, na.rm = TRUE)
    IPRESscore <- colMeans(z_score)

    return(returnAsInput(userdata = dataset, result = IPRESscore,
                         SignName = "IPRES", datasetm))
}

#'  CIS (carcinoma-in situ) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the CIS score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
CISSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "bladder"){

    firstCheck(nametype, tumorTissue, "CISSign")

    if(nametype!="SYMBOL"){
        CISdata$Gene_Symbol  <- mapIds(
            org.Hs.eg.db, keys = CISdata$Gene_Symbol, column = nametype,
            keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    Signature_up <- CISdata[grep('CISup', CISdata$Category),]
    Signature_down <- CISdata[grep('CISdown', CISdata$Category),]

    percentageOfGenesUsed("CISSign", datasetm, Signature_up$Gene_Symbol, "up")
    percentageOfGenesUsed("CISSign", datasetm, Signature_down$Gene_Symbol, "down")

    med_data_up <- colMeans(log2(datasetm[intersect(
        row.names(datasetm), Signature_up$Gene_Symbol),]))
    med_data_down <- colMeans(log2(datasetm[intersect(
        row.names(datasetm), Signature_down$Gene_Symbol),]))

    CISscore <- med_data_up - med_data_down

    return(returnAsInput(
        userdata = dataset, result = CISscore, SignName = "CIS", datasetm))
}

#' Glycolysis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the glycolysis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
glycolysisSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "breast",
                           author = "Jiang"){

    firstCheck(nametype, tumorTissue, "glycolysisSign", author)

    Glycolysisdata <- get(paste0("Glycolysis", author, "data"))
    datasetm <- getMatrix(dataset)

    score <- coefficientsScore(Glycolysisdata, datasetm,
                                nametype, "glycolysisSign")

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = paste0("Glycolysis", author), datasetm))
}

#' Autophagy Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Autophagy score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
autophagySign <- function(dataset, nametype = "SYMBOL",
                          tumorTissue = "brain", author = "Xu"){

    firstCheck(nametype, tumorTissue, "autophagySign", author)

    Autophagydata <- get(paste0("Autophagy", author, "data"))
    datasetm <- getMatrix(dataset)
    score <- coefficientsScore(Autophagydata, datasetm,
                                nametype, "autophagySign")

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = paste0("Autophagy", author), datasetm))
}

#' Extracellular Matrix Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the ECM scores
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ECMSign <- function(dataset, nametype = "SYMBOL",
                    tumorTissue = "pan-tissue", pvalues = FALSE,
                    nperm = 100, ...){

    firstCheck(nametype, tumorTissue, "ECMSign")

    if(nametype!="SYMBOL"){
        ECMdata$Gene_Symbol <- mapIds(
            org.Hs.eg.db, keys = ECMdata$Gene_Symbol,
            column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    Signature_up <- ECMdata[grep('ECMup', ECMdata$Category),]
    Signature_down <- ECMdata[grep('ECMdown', ECMdata$Category),]

    percentageOfGenesUsed("ECMSign", datasetm, Signature_up$Gene_Symbol, "up")
    percentageOfGenesUsed("ECMSign", datasetm, Signature_down$Gene_Symbol, "down")

    gene_sets <- list(ECMup = Signature_up$Gene_Symbol,
                      ECMdown = Signature_down$Gene_Symbol)

    dots <- list(...)

    args <- matchArguments(dots,list(
        expr = datasetm, gset.idx.list = gene_sets,
        method = "gsva", kcdf = "Gaussian", min.sz = 5,
        ssgsea.norm = FALSE, verbose = FALSE))

    gsva_count <- suppressWarnings(do.call(gsva, args))

    if(pvalues){
        gsva_pval <- GSVAPvalues(
            expr = datasetm, gset.idx.list = gene_sets,
            gsvaResult = gsva_matrix, nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(
        userdata = dataset, result = gsva_count, SignName = "", datasetm))
}


#' Homologous Recombination Deficiency Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the HRDS scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats setNames
#'
#' @import org.Hs.eg.db
#'
#' @export
HRDSSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "HRDSSign")

    if(nametype!="SYMBOL"){
        HRDSdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys= HRDSdata$Gene_Symbol,
                                       column = nametype, keytype = "SYMBOL",
                                       multiVals = "first")}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("HRDSSign", datasetm, HRDSdata$Gene_Symbol)

    HRDS_P <- datasetm[intersect(row.names(datasetm), HRDSdata[
            HRDSdata$correlation == 1, ]$Gene_Symbol), ]
    HRDS_N <- datasetm[intersect(row.names(datasetm), HRDSdata[
            HRDSdata$correlation == -1, ]$Gene_Symbol), ]

    score <- unlist(lapply(seq_len(ncol(datasetm)), function(x){
        tmp <- t.test(HRDS_P[,x], HRDS_N[,x], alternative = "two.sided")
        tmp[["statistic"]]}))

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = "HRDS", datasetm))
}


#' Adult Intestinal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the ISC scores will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ISCSign <- function(dataset, nametype= "SYMBOL", tumorTissue = "colon"){

    firstCheck(nametype, tumorTissue, "ISCSign")

    if(nametype!= "SYMBOL"){
        ISCdata <- lapply(ISCdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype,
                                    keytype="SYMBOL", multiVals="first")))}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ISCSign", datasetm, ISCdata$Lgr5, "Lgr5_ISC")
    percentageOfGenesUsed("ISCSign", datasetm, ISCdata$Eph, "EphB2_ISC")
    percentageOfGenesUsed("ISCSign", datasetm, ISCdata$TA, "Late_TA")
    percentageOfGenesUsed("ISCSign", datasetm, ISCdata$pro, "Proliferation")

    ISCscores <- sapply(ISCdata, function(x)
        colMeans(datasetm[intersect(
            x, row.names(datasetm - colMeans(datasetm))),]))

    return(returnAsInput(userdata = dataset, result = t(ISCscores),
                         SignName = "ISC", datasetm))
}

#' VEGF Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the VEGF score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
VEGFSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "VEGFSign")

    datasetm <- getMatrix(dataset)

    score <- log2(statScore(VEGFdata, datasetm, nametype, "mean", "VEGFSign"))

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = "VEGF", datasetm))
}

#' Angiogenesis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Angiogenesis score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
angioSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "angioSign")

    datasetm <- getMatrix(dataset)

    score <-statScore(
        Angiogenesisdata, datasetm, nametype, typeofstat = "median",
        namesignature = "angioSign")

    return(returnAsInput(userdata = dataset, result = as.vector(scale(score)),
                         SignName = "Angiogenesis", datasetm))
}

#' DNA Repair Signature
#'
##' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Angiogenesis score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
DNArepSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "DNArepSign")

    if(nametype!="SYMBOL"){DNArepairdata$DNAre <- mapIds(
        org.Hs.eg.db, keys = DNArepairdata$DNAre, column = nametype,
        keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("DNArepSign", datasetm, DNArepairdata$DNAre)

    DNAdatahigh <- DNArepairdata[DNArepairdata$status=="high",]
    DNAdatalow <- DNArepairdata[DNArepairdata$status=="low",]

    medianexp <- apply(datasetm, 2, median)

    DNArepscore <- apply(
        datasetm[DNAdatahigh$DNAre, ] > medianexp, 2, sum,na.rm = TRUE)+
        apply(datasetm[DNAdatalow$DNAre, ] > medianexp, 2, sum, na.rm = TRUE)

    return(returnAsInput(userdata = dataset, result = DNArepscore,
                         SignName = "DNArepair", datasetm))
}
