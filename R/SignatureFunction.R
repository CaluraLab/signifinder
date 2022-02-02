
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
#' @param author first author of the specific signature publication.
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the Epithelial and Mesenchymal
#' scores are added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#'
#' @export
EMTSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray",
    author = "Miow", pvalues = FALSE, nperm = 100, ...) {

    consistencyCheck(nametype, "EMTSign", author)

    datasetm <- getMatrix(dataset)

    if(author == "Miow"){

        EMT_Miow$SYMBOL <- geneIDtrans(nametype, EMT_Miow$SYMBOL)

        EL <- EMT_Miow[grep('Epithelial-like', EMT_Miow$class),]
        ML <- EMT_Miow[grep('Mesenchymal-like', EMT_Miow$class),]

        percentageOfGenesUsed("EMTSign", datasetm, EL$SYMBOL, "epithelial")
        percentageOfGenesUsed("EMTSign", datasetm, ML$SYMBOL, "mesenchymal")

        gene_sets <- list(Epithelial = EL$SYMBOL, Mesenchymal = ML$SYMBOL)

        dots <- list(...)
        kcdftype <- ifelse(inputType == "microarray", "Gaussian", "Poisson")
        args <- matchArguments(dots, list(
            expr = datasetm, gset.idx.list = gene_sets, method = "ssgsea",
            kcdf = kcdftype, min.sz = 5, ssgsea.norm = FALSE, verbose = FALSE))
        gsva_matrix <- suppressWarnings(do.call(gsva, args))

        if(pvalues){
            gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets,
                gsvaResult = gsva_matrix, nperm = nperm, args = args)
            gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

        return(returnAsInput(
            userdata = dataset, result = gsva_matrix, SignName = "", datasetm))

    } else {
        if(author == "Mak") {

            EMT_Mak$SYMBOL <- geneIDtrans(nametype, EMT_Mak$SYMBOL)

            Sign_E <- EMT_Mak$SYMBOL[EMT_Mak$class=="E"]
            Sign_M <- EMT_Mak$SYMBOL[EMT_Mak$class=="M"]

            percentageOfGenesUsed("EMTSign", datasetm, Sign_E, "epithelial")
            percentageOfGenesUsed("EMTSign", datasetm, Sign_M, "mesenchymal")

            columnNA <- managena(datasetm, genes = EMT_Mak$SYMBOL)
            score <- colMeans(
                datasetm[intersect(Sign_M, row.names(datasetm)),]
                ) - colMeans(
                datasetm[intersect(Sign_E, row.names(datasetm)),])
            score[columnNA > 0.9] <- NA

        } else if(author == "Cheng") {

            EMT_Cheng$SYMBOL <- geneIDtrans(nametype, EMT_Cheng$SYMBOL)

            percentageOfGenesUsed("EMTSign", datasetm, EMT_Cheng$SYMBOL)

            datasetm_n <- datasetm[intersect(
                row.names(datasetm), EMT_Cheng$SYMBOL), ]
            datasetm_n <- if(inputType == "rnaseq") {log2(datasetm_n)
                } else {datasetm_n}
            columnNA <- managena(datasetm = datasetm_n,
                                genes = EMT_Cheng$SYMBOL)
            score <- prcomp(t(datasetm_n))$x[,1]
            score[columnNA > 0.9] <- NA
        }
        return(returnAsInput(userdata = dataset, result = score,
                             SignName = paste0("EMT_", author), datasetm))}
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
#' @export
pyroptosisSign <- function(dataset, nametype = "SYMBOL", inputType = "rnaseq",
                author = "Ye", hgReference = "hg19"){

    consistencyCheck(nametype, "pyroptosisSign", author)

    Pyroptosisdata <- get(paste0("Pyroptosis_", author))

    datasetm <- getMatrix(dataset)

    if(author == "Ye"){
        # datasetm_n <- scale(datasetm)
        datasetm_n <- dataTransformation(datasetm, "FPKM", hgReference)
    } else if (author == "Shao"){
        if(inputType == "rnaseq"){
            datasetm_n <- dataTransformation(datasetm, "FPKM", hgReference)
        } else {datasetm_n <- datasetm}
    } else if (author == "Lin"){
        datasetm_n <- dataTransformation(datasetm, "TPM", hgReference)
    } else {datasetm_n <- datasetm}

    score <- coeffScore(Pyroptosisdata, datasetm_n, nametype, "pyroptosisSign")

    return(returnAsInput(
        userdata = dataset, result = score,
        SignName = paste0("Pyroptosis_", author), datasetm))
}


#' Ferroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the Ferroptosis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @export
ferroptosisSign <- function(dataset, nametype = "SYMBOL", inputType = "rnaseq",
                    author = "Ye", hgReference = "hg19"){

    consistencyCheck(nametype, "ferroptosisSign", author)
    Ferrdata <- get(paste0("Ferroptosis_", author))

    datasetm <- getMatrix(dataset)

    if(author == "Liu"){
        if(inputType == "rnaseq"){
            datasetm_n <- dataTransformation(datasetm, "FPKM", hgReference)
        } else {datasetm_n <- datasetm}
    } else if (author == "Li"){
        datasetm_n <- datasetm
    } else {
        datasetm_n <- datasetm}

    score <- coeffScore(Ferrdata, datasetm_n, nametype, "ferroptosisSign")

    if(author == "Liang" ){score <- exp(score)}

    return(returnAsInput(userdata = dataset, result = score,
        SignName = paste0("Ferroptosis_", author), datasetm))
}


#' Lipid Metabolism Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Lipid scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @export
lipidMetabolismSign <- function(dataset, nametype = "SYMBOL") {

    consistencyCheck(nametype, "lipidMetabolismSign")

    datasetm <- getMatrix(dataset)
    score <- coeffScore(LipidMetabolism_Zheng, datasetm,
                        nametype, "lipidMetabolismSign")

    return(returnAsInput(userdata = dataset, result = score,
                SignName = "LipidMetabolism_Zheng", datasetm))
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
#' @import org.Hs.eg.db
#'
#' @export
hypoxiaSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray"){

    consistencyCheck(nametype, "hypoxiaSign")

    # if(nametype=="ENSEMBL") {genetouse <- Hypoxiadata$Gene_Ensembl}

    datasetm <- getMatrix(dataset)
    datasetm_n <- if(inputType == "rnaseq") {log2(datasetm)} else {datasetm}
    score <- statScore(Hypoxia_Buffa$SYMBOL, datasetm = datasetm_n,
                        nametype = "SYMBOL", typeofstat = "median",
                        namesignature = "hypoxiaSign")

    return(returnAsInput(userdata = dataset, result = as.vector(scale(score)),
        SignName = "Hypoxia_Buffa", datasetm))
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
#'
#' @export
platinumResistanceSign <- function(dataset, nametype = "SYMBOL",
        pvalues = FALSE, nperm = 100, ...){

    consistencyCheck(nametype, "platinumResistanceSign")

    PlatinumResistancedata$Gene_Symbol <- geneIDtrans(
        nametype, PlatinumResistancedata$Gene_Symbol)

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
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the Immunogenic scores will
#' be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom labstatR meang
#'
#' @export
immunoScoreSign <- function(dataset, nametype = "SYMBOL", author = "Hao",
                            inputType = "rnaseq", hgReference = "hg19"){

    consistencyCheck(nametype, "immunoScoreSign", author)

    datasetm <- getMatrix(dataset)

    if(author == "Hao"){

        ImmunoScoreHaodata$genes <- geneIDtrans(
            nametype, ImmunoScoreHaodata$genes)

        g <- intersect(row.names(datasetm), ImmunoScoreHaodata$genes)

        percentageOfGenesUsed("immunoScoreSign", datasetm,
                              ImmunoScoreHaodata$genes)

        datasetm_n <- datasetm[g,]
        if(inputType == "rnaseq"){
            datasetm_n <- log2(dataTransformation(
                datasetm_n, "FPKM", hgReference) + 0.01)}

        ImmunoScoreHaodata <- ImmunoScoreHaodata[ImmunoScoreHaodata$genes%in%g,]
        columnNA <- managena(datasetm_n, g)
        SE <- (ImmunoScoreHaodata$HR - ImmunoScoreHaodata$`95CI_L`)/1.96
        k <- (1 - ImmunoScoreHaodata$HR)/SE
        score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(p){
                sum(k*datasetm_n[,p], na.rm = TRUE)}))
        score[columnNA > 0.9] <- NA
    } else if(author == "Roh"){
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
consensusOVSign <- function(dataset, nametype = "SYMBOL",
                            method = "consensusOV", ...){

    consistencyCheck(nametype, "consensusOVSign")

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
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the
#' IPS, MHC, CP, EC and SC scores will be added in the
#' \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom stats sd
#'
#' @export
IPSSign <- function(dataset, nametype = "SYMBOL", hgReference = "hg19"){

    consistencyCheck(nametype, "IPSSign")

    IPSdata[,c(1,2)] <- data.frame(lapply(IPSdata[,c(1,2)], function(x){
        geneIDtrans(nametype, x)}))
    colnames(IPSdata) <- c("GENE","NAME","CLASS","WEIGHT")

    datasetm <- getMatrix(dataset)
    sample_names <- colnames(datasetm)

    percentageOfGenesUsed("IPSSign", datasetm, IPSdata$GENE)

    ind <- which(is.na(match(IPSdata$GENE, row.names(datasetm))))
    MISSING_GENES <- IPSdata$GENE[ind]
    if (length(MISSING_GENES)>0) {cat("Differently named or missing genes: ",
                                    MISSING_GENES, "\n")}

    datasetm_n <- log2(dataTransformation(datasetm, "TPM", hgReference) + 1)

    IPS <- NULL; MHC <- NULL; CP <- NULL; EC <- NULL; SC <- NULL; AZ <- NULL
    for (i in 1:length(sample_names)) {
        GE <- datasetm_n[,i]
        Z1 <- (datasetm_n[match(
            IPSdata$GENE, row.names(datasetm_n)),i]-mean(GE))/sd(GE)
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
#' @export
matrisomeSign <- function(dataset, nametype = "SYMBOL") {

    consistencyCheck(nametype, "matrisomeSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(Matrisome_Yuzhalin$SYMBOL, datasetm, nametype,
                        "median", "matrisomeSign")

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "Matrisome_Yuzhalin", datasetm))
}


#' Mitotic Index
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @export
mitoticIndexSign <- function(dataset, nametype = "SYMBOL") {

    consistencyCheck(nametype, "mitoticIndexSign")

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
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom labstatR meang
#'
#' @export
ImmuneCytSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray",
            author = "Rooney", hgReference = "hg19"){

    consistencyCheck(nametype, "ImmuneCytSign", author)

    datasetm <- getMatrix(dataset)
    if(author == "Rooney"){
        if(inputType == "rnaseq"){
            datasetm_n <- dataTransformation(datasetm, "TPM", hgReference)+0.01
        } else {datasetm_n <- datasetm}
         score <- statScore(
             ImmuneCytRooneydata, datasetm_n, nametype, "meang", "ImmuneCytSign")
    } else if(author == "Davoli") {
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
#' @export
IFNSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "IFNSign")

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
#' @export
expandedImmuneSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "expandedImmuneSign")

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
#' @importFrom labstatR meang
#'
#' @export
TLSSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "TLSSign")

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
#' @export
CD49BSCSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "CD49BSCSign")

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
#' @export
CINSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "CINSign")

    datasetm <- getMatrix(dataset)

    score25 <- statScore(CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"],
                        scale(datasetm, center = TRUE, scale = FALSE),
                        nametype, "sum", "CINSign")
    score70 <- statScore(CIN_Carter$SYMBOL,
                        scale(datasetm, center = TRUE, scale = FALSE),
                        nametype, "sum", "CINSign")
    score <- data.frame(CIN25 = score25, CIN70 = score70)

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
#' @export
CCSSign <- function(dataset, nametype = "SYMBOL", author = "Lundberg"){

    consistencyCheck(nametype, "CCSSign", author)

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
#' @export
chemokineSign <- function(dataset, nametype = "SYMBOL",
                          inputType = "microarray"){

    consistencyCheck(nametype, "chemokineSign")

    Chemokines_Messina$SYMBOL <- geneIDtrans(nametype, Chemokines_Messina$SYMBOL)

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ChemokineSign", datasetm, Chemokines_Messina$SYMBOL)

    datasetm_n <- datasetm[intersect(row.names(datasetm), Chemokines_Messina$SYMBOL), ]
    datasetm_n <- if(inputType == "rnaseq") {log2(datasetm_n)
    } else {datasetm_n}
    columnNA <- managena(datasetm_n, Chemokines_Messina$SYMBOL)
    score <- prcomp(t(datasetm), center = TRUE, scale = TRUE)$x[, 1]
    score[columnNA > 0.9] <- NA

    return(returnAsInput(userdata = dataset, result = score,
                        SignName = "Chemokines_Messina", datasetm))
}

#' Adult Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @export
ASCSign <- function(dataset, nametype= "SYMBOL"){

    consistencyCheck(nametype, "ASCSign")

    ASCdata <- geneIDtrans(nametype, ASCdata)

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ASCSign", datasetm, ASCdata)

    datasetm_n <- log2(datasetm[row.names(datasetm) %in% ASCdata,] + 1)

    columnNA <- managena(datasetm_n, ASCdata)

    score <- rowSums(
        scale(t(datasetm_n), center = TRUE, scale = TRUE), na.rm = TRUE)
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
#'
#' @export
PassONSign <- function(dataset, nametype = "SYMBOL", ...){

    consistencyCheck(nametype, "PassONSign")

    PASS.ONdata <- lapply(PASS.ONdata, function(x){geneIDtrans(nametype, x)})

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
#'
#' @export
IPRESSign <- function(dataset, nametype = "SYMBOL",
                      pvalues = FALSE, nperm = 100, ...) {

    consistencyCheck(nametype, "IPRESSign")

    IPRESdata <- lapply(IPRESdata, function(x){geneIDtrans(nametype, x)})

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
#' @export
CISSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "CISSign")

    CISdata$Gene_Symbol <- geneIDtrans(nametype, CISdata$Gene_Symbol)

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
#' @export
glycolysisSign <- function(dataset, nametype = "SYMBOL", author = "Jiang"){

    consistencyCheck(nametype, "glycolysisSign", author)

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
#' @export
autophagySign <- function(dataset, nametype = "SYMBOL", author = "Xu"){

    consistencyCheck(nametype, "autophagySign", author)

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
#'
#' @export
ECMSign <- function(dataset, nametype = "SYMBOL",
                    pvalues = FALSE, nperm = 100, ...){

    consistencyCheck(nametype, "ECMSign")

    ECMdata$Gene_Symbol <- geneIDtrans(nametype, ECMdata$Gene_Symbol)

    datasetm <- getMatrix(dataset)

    Signature_up <- ECMdata[grep('ECMup', ECMdata$Category),]
    Signature_down <- ECMdata[grep('ECMdown', ECMdata$Category),]

    percentageOfGenesUsed("ECMSign", datasetm, Signature_up$Gene_Symbol, "up")
    percentageOfGenesUsed("ECMSign", datasetm, Signature_down$Gene_Symbol, "down")

    gene_sets <- list(ECMup = Signature_up$Gene_Symbol,
                      ECMdown = Signature_down$Gene_Symbol)

    dots <- list(...)

    args <- matchArguments(dots,list(expr = datasetm,
                                     gset.idx.list = gene_sets,
                                     method = "ssgsea", kcdf = "Poisson",
                                     min.sz = 5, ssgsea.norm = FALSE,
                                     verbose = FALSE))

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
#' @export
HRDSSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "HRDSSign")

    HRDSdata$Gene_Symbol <- geneIDtrans(nametype, HRDSdata$Gene_Symbol)

    datasetm <- getMatrix(dataset)
    datasetm_n <- datasetm - rowMedians(datasetm)

    percentageOfGenesUsed("HRDSSign", datasetm, HRDSdata$Gene_Symbol)

    HRDS_P <- datasetm_n[intersect(row.names(datasetm_n), HRDSdata[
            HRDSdata$correlation == 1, ]$Gene_Symbol), ]
    HRDS_N <- datasetm_n[intersect(row.names(datasetm_n), HRDSdata[
            HRDSdata$correlation == -1, ]$Gene_Symbol), ]

    score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(x){
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
#' @export
ISCSign <- function(dataset, nametype= "SYMBOL"){

    consistencyCheck(nametype, "ISCSign")

    ISCdata <- lapply(ISCdata, function(x){geneIDtrans(nametype, x)})

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
#' @export
VEGFSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "VEGFSign")

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm)

    score <- statScore(VEGFdata, datasetm, nametype, "mean", "VEGFSign")

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
#' @export
angioSign <- function(dataset, nametype = "SYMBOL"){

    consistencyCheck(nametype, "angioSign")

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
#' @export
DNArepSign <- function(dataset, nametype = "SYMBOL", inputType = "microarray"){

    consistencyCheck(nametype, "DNArepSign")

    DNArepairdata$DNAre <- geneIDtrans(nametype, DNArepairdata$DNAre)

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("DNArepSign", datasetm, DNArepairdata$DNAre)

    datasetm_n <- if(inputType == "rnaseq") {log2(datasetm)
        } else {datasetm}
    datasetm_n <- datasetm_n[row.names(datasetm_n) %in% DNArepairdata$DNAre, ]
    datasetm_n <- scale(t(datasetm_n), center = TRUE, scale = FALSE)

    medianexp <- apply(datasetm_n, 2, median)

    genes_h <- intersect(
    colnames(datasetm_n), DNArepairdata[DNArepairdata$status=="high",]$DNAre)
    genes_l <- intersect(
    colnames(datasetm_n), DNArepairdata[DNArepairdata$status=="low",]$DNAre)

    score <- rowSums(datasetm_n[,genes_h] > medianexp[genes_h]) +
        rowSums(datasetm_n[,genes_l] < medianexp[genes_l])

    return(returnAsInput(userdata = dataset, result = score,
                         SignName = "DNArepair", datasetm))
}
