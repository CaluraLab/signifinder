
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

    eper <- (sum(Signature_EL$Gene_Symbol %in% row.names(datasetm))/
                nrow(Signature_EL))*100
    mper <- (sum(Signature_ML$Gene_Symbol %in% row.names(datasetm))/
                nrow(Signature_ML))*100
    cat(paste0("EMTSign function is using ", round(eper),
                "% of epithelial-like genes\n", "EMTSign function is using ",
                round(mper), "% of mesenchymal-like genes\n"))

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

    } else if(tumorTissue=="pan-tissue" & author == "Mak"){
        if(nametype!="SYMBOL"){
            EMTMakdata$Gene_Symbol <- mapIds(
                org.Hs.eg.db, keys = EMTMakdata$Gene_Symbol, column = nametype,
                keytype = "SYMBOL", multiVals = "first")}

        Sign_E <- EMTMakdata[grep('E', EMTMakdata$Category),]
        Sign_M <- EMTMakdata[grep('M', EMTMakdata$Category),]

        eper <- (sum(Sign_E$Gene_Symbol %in% row.names(datasetm))/
                     nrow(Sign_E))*100
        mper <- (sum(Sign_M$Gene_Symbol %in% row.names(datasetm))/
                     nrow(Sign_M))*100
        cat(paste0("EMTSign function is using ", round(eper),
                   "% of epithelial-like genes\n",
                   "EMTSign function is using ", round(mper),
                   "% of mesenchymal-like genes\n"))

        if(inputType == "rnaseq"){datasetm <- log2(datasetm)}
        EMTscore <- colMeans(datasetm[intersect(
            Sign_M$Gene_Symbol, row.names(datasetm)),]) - colMeans(datasetm[
                intersect(Sign_E$Gene_Symbol, row.names(datasetm)), ])

    } else if(tumorTissue=="breast" & author == "Cheng"){

        if(nametype!="SYMBOL") {
            EMTChengdata <- mapIds(org.Hs.eg.db, keys = EMTChengdata,
                                   column = nametype, keytype = "SYMBOL",
                                   multiVals = "first")}

        Ebper <- (sum(EMTChengdata %in% row.names(datasetm))/
                      length(EMTChengdata))*100
        cat(paste0("EMTSign function is using ", round(Ebper),
                   "% of genes\n"))

        datasetm <- datasetm[intersect(row.names(datasetm), EMTChengdata), ]
        if(inputType == "rnaseq"){datasetm <- log2(datasetm)}
        # managena(datasetm, EMTChengdata)
        EMTscore <- prcomp(t(datasetm))$x[,1]}

    return(returnAsInput(userdata = dataset, result = EMTscore,
                         SignName = paste0("EMT", author), datasetm))
}


#' Pyroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Pyroptosis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
pyroptosisSign <- function(dataset, nametype = "SYMBOL", inputType = "rnaseq",
                            tumorTissue = "ovary", author = "Ye"){

    firstCheck(nametype, tumorTissue, "pyroptosisSign", author)

    Pyroptosisdata <- get(paste0("Pyroptosis", author, "data"))

    datasetm <- getMatrix(dataset)

    if(tumorTissue == "ovary" & author == "Ye"){
        datasetm <- scale(datasetm)
        datasetm <- dataTransformation(datasetm, "FPKM")
    } else if (tumorTissue == "stomach" & author == "Shao"){
        if(inputType == "rnaseq"){
            datasetm <- dataTransformation(datasetm, "FPKM")}
    } else if (tumorTissue == "lung" & author == "Lin"){
        datasetm <- dataTransformation(datasetm, "TPM")
    }

    Piroscore <- coefficientsScore(
        Pyroptosisdata, datasetm = datasetm, nametype = nametype,
        namesignature = "pyroptosisSign")

    return(returnAsInput(
        userdata = dataset, result = Piroscore,
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
    ferrscore <- coefficientsScore(Ferroptosisdata, datasetm = datasetm,
                                    nametype = nametype,
                                    namesignature = "ferroptosisSign")

    if(tumorTissue == "liver" & author == "Liang" ){ferrscore <- exp(ferrscore)}

    return(returnAsInput(userdata = dataset, result = ferrscore,
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
    lipidscore <- coefficientsScore(LipidMetabolismdata, datasetm = datasetm,
                                    nametype = nametype,
                                    namesignature = "lipidMetabolismSign")

    return(returnAsInput(userdata = dataset, result = lipidscore,
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
    if(inputType == "rnaseq"){datasetm <- log2(datasetm)}
    score <- statScore(genetouse, datasetm = datasetm, nametype = "SYMBOL",
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
                                    tumorTissue = "ovary", pvalues = FALSE,
                                    nperm = 100, ...){

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

    upper <- (sum(Signature_up$Gene_Symbol %in% row.names(datasetm))/
                 nrow(Signature_up))*100
    downper <- (sum(Signature_down$Gene_Symbol %in% row.names(datasetm))/
                 nrow(Signature_down))*100
    cat(paste0("platinumResistanceSign function is using ", round(upper),
                "% of up-genes\n", "platinumResistanceSign function is using ",
                round(downper), "% of down-genes\n"))

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

    if(tumorTissue=="ovary"){
        if(nametype!="SYMBOL"){
            ImmunoScoreHaodata$genes <- mapIds(
                org.Hs.eg.db, keys = ImmunoScoreHaodata$genes,
                column = nametype, keytype = "SYMBOL", multiVals = "first")}

        g <- intersect(row.names(datasetm), ImmunoScoreHaodata$genes)
        gper <- (length(g)/length(ImmunoScoreHaodata$genes))*100
        cat(paste0(
            "immunoScoreSign function is using ", round(gper), "% of genes\n"))

        subdataset <- datasetm[g,]
        ImmunoScoreHaodata <- ImmunoScoreHaodata[ImmunoScoreHaodata$genes%in%g,]

        SE <- (ImmunoScoreHaodata$HR - ImmunoScoreHaodata$`95CI_L`)/1.96
        k <- (1 - ImmunoScoreHaodata$HR)/SE
        ImmunoScores <- unlist(lapply(
            seq_len(ncol(subdataset)), function(p){
                sum(k*subdataset[,p], na.rm = TRUE)}))
    } else if(tumorTissue=="pan-tissue"){
        ImmunoScores <- statScore(
            ImmunoScoreRohdata, datasetm = datasetm, nametype = nametype,
            typeofstat = "meang", namesignature = "immunoScoreSign")}

    return(returnAsInput(
        userdata = dataset, result = ImmunoScores,
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
        datasetm <- datasetm[!is.na(genename),]
        genename <- genename[!is.na(genename)]
        datasetm <- datasetm[!duplicated(genename),]
        genename <- genename[!duplicated(genename)]
    } else {genename <- row.names(datasetm)}

    consensus_subtypes <- get.subtypes(expression.dataset=datasetm,
                                    entrez.ids=genename, method=method, ...)

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

    ipsper <- (sum(IPSdata$GENE %in% row.names(datasetm))/nrow(IPSdata))*100
    cat(paste0("IPSSign function is using", round(ipsper), "% of genes\n"))

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
    median_cm <- statScore(Matrisomedata, datasetm = datasetm,
                        nametype = nametype, typeofstat = "median",
                        namesignature = "matrisomeSign")

    return(returnAsInput(userdata = dataset, result = median_cm,
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
    MI_means <- statScore(MitoticIndexdata, datasetm = datasetm,
                        nametype = nametype, typeofstat = "mean",
                        namesignature = "mitoticIndexSign")

    return(returnAsInput(userdata = dataset, result = MI_means,
                        SignName = "MitoticIndex", datasetm))
}


#' Local Immune Cytolytic Activity (CYT) Signature
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
CYTSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "CYTSign")

    datasetm <- getMatrix(dataset)
    CYTScore <- statScore(CYTdata, datasetm = datasetm, nametype = nametype,
                        typeofstat = "meang", namesignature = "CYTSign")

    return(returnAsInput(userdata = dataset, result = CYTScore,
                        SignName = "CYT", datasetm))
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
    IFNscore <- statScore(IFNdata, datasetm = datasetm, nametype = nametype,
                        typeofstat = "mean", namesignature = "IFNSign")

    return(returnAsInput(userdata = dataset, result = IFNscore,
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
    ExpandedImmunescore <- statScore(ExpandedImmunedata, datasetm = datasetm,
                                    nametype = nametype, typeofstat = "mean",
                                    namesignature = "expandedImmuneSign")

    return(returnAsInput(userdata = dataset, result = ExpandedImmunescore,
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
    TLSScore <- statScore(TLSdata, datasetm = datasetm, nametype = nametype,
                        typeofstat = "meang", namesignature = "TLSSign")

    return(returnAsInput(userdata = dataset, result = TLSScore,
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
    CD49BSCScore <- coefficientsScore(CD49BSCdata, datasetm = datasetm,
                    nametype = nametype, namesignature = "CD49BSCSign")

    return(returnAsInput(userdata = dataset, result = CD49BSCScore,
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
                    tumorTissue = "breast, lung, brain, lymphatic system",
                    typeofCIN = 70){

    firstCheck(nametype, tumorTissue, "CINSign")
    if (!(typeofCIN %in% c(70,25))){stop("typeofCIN must be either 70 or 25")}

    CINdata <- CINdata[seq_len(typeofCIN)]

    datasetm <- getMatrix(dataset)

    CINscore <- statScore(
        CINdata, datasetm = datasetm, nametype = nametype, typeofstat = "sum",
        namesignature = "CINSign")

    return(returnAsInput(userdata = dataset, result = CINscore,
                        SignName = "CIN", datasetm))
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

    if(tumorTissue=="pan-tissue"&author=="Lundberg"){

        CCSscore <- statScore(
            CCSLundbergdata, datasetm = datasetm, nametype = nametype,
            typeofstat = "sum", namesignature = "CCSSign")

    } else if(tumorTissue=="pan-tissue"&author == "Davoli"){

        CCSscore <- statScore(
            CCSDavolidata, datasetm = datasetm, nametype = nametype,
            typeofstat = "mean", namesignature = "CCSSign")}

    return(returnAsInput(
        userdata = dataset, result = CCSscore, SignName = paste0("CCS", author),
        datasetm))
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

    if(nametype!="SYMBOL") {Chemokinedata <- mapIds(org.Hs.eg.db,
                                                    keys = Chemokinedata,
                                                    column = nametype,
                                                    keytype = "SYMBOL",
                                                    multiVals = "first")}

    datasetm <- getMatrix(dataset)

    chemoper <- (sum(Chemokinedata %in% row.names(datasetm))/
                     length(Chemokinedata))*100
    cat(paste0("ChemokineSign function is using ",
               round(chemoper), "% of genes\n"))

    datasetm <- datasetm[intersect(row.names(datasetm), Chemokinedata), ]
    z_score <- (datasetm - rowMeans(datasetm))/
        sapply(as.data.frame(t(datasetm)), sd, na.rm = TRUE)
    pc1_score <- prcomp(t(z_score))$x[,1]
    return(returnAsInput(userdata = dataset, result = pc1_score,
                                       SignName = "Chemokines", datasetm))
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
                              tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "ASCSign")

    if(nametype!="SYMBOL"){
        ASCdata <- mapIds(org.Hs.eg.db, keys = ASCdata, column = nametype,
                          keytype = "SYMBOL", multiVals = "first")}
    datasetm <- getMatrix(dataset)
    ASCper <- (sum(ASCdata %in% row.names(datasetm))/length(ASCdata))*100
    cat(paste0("AdultStemCellSign function is using ",
               round(ASCper), "% of AdultStemCell genes\n"))

    ASCdatasetm <- log2(datasetm[row.names(datasetm) %in% ASCdata, ] + 1)
    ASCscore <- colSums((ASCdatasetm - rowMeans(ASCdatasetm))/
                            sapply(as.data.frame(t(ASCdatasetm)), sd,
                                   na.rm = TRUE))

    return(returnAsInput(userdata = dataset, result = ASCscore,
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

    passper <- (sum(unlist(PASS.ONdata) %in% row.names(datasetm))/
                    sum(lengths(PASS.ONdata)))*100

    cat(paste0("passONsignature function is using ", round(passper),
               "% of passON-related genes\n"))

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm,
                                      gset.idx.list = PASS.ONdata,
                                      method = "ssgsea",
                                      kcdf = "Poisson",
                                      min.sz = 5,
                                      ssgsea.norm = TRUE,
                                      verbose = TRUE))

    gsva_matrix <- suppressWarnings(do.call(gsva, args))

    gsva_mean <- sapply(X = colnames(gsva_matrix) , FUN = function(x) {
        weighted.mean(gsva_matrix[,x], c(23, 26, 81, 18))
    })


    return(returnAsInput(userdata = dataset, result = gsva_mean,
                                       SignName = "PassON", datasetm))
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
        ipresdata <- mapIds(org.Hs.eg.db, keys = ipresdata,
                            column = nametype, keytype = "SYMBOL",
                            multiVals = "first")}

    datasetm <- getMatrix(dataset)


    passper <- (sum(unlist(ipresdata) %in% row.names(datasetm))/
                    sum(lengths(ipresdata)))*100

    cat(paste0("IPRESsignature function is using ", round(passper),
               "% of IPRES-related genes\n"))


    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm,
                 gset.idx.list = ipresdata,
                 method = "ssgsea", kcdf = "Poisson",
                 min.sz = 5, ssgsea.norm = FALSE,
                 verbose = FALSE))

    gsva_matrix <- suppressWarnings(do.call(gsva, args))

    if(pvalues){
        gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets,
                                 gsvaResult = gsva_matrix,
                                 nperm = nperm, args = args)
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
        CISdata  <- mapIds(org.Hs.eg.db, keys = CISdata, column = nametype,
                           keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    upper <- (sum(CISdata$UP %in% row.names(datasetm))/
                  length(CISdata$UP))*100
    downper <- (sum(CISdata$DOWN %in% row.names(datasetm))/
                    length(CISdata$DOWN))*100
    cat(paste0("CISSign function is using ", round(upper),
               "% of up-genes\n", "CISSign function is using ",
               round(downper), "% of down-genes\n"))

    med_data_up <- colMeans(log2(datasetm[intersect(
        row.names(datasetm), CISdata$UP),]))
    med_data_down <- colMeans(log2(datasetm[intersect(
        row.names(datasetm), CISdata$DOWN),]))

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

    Glycolisisdata <- get(paste0("Glycolisis", author, "data"))

    datasetm <- getMatrix(dataset)

    Glycoscore <- coefficientsScore(Glycolisisdata, datasetm = datasetm,
                                   nametype = nametype,
                                   namesignature = "glycolysisSign")

    return(returnAsInput(userdata = dataset, result = Glycoscore,
                         SignName = paste0("Glycolysis",tumorTissue), datasetm))
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

    Autoscore <- coefficientsScore(Autophagydata, datasetm = datasetm,
                 nametype = nametype, namesignature = "autophagySign")

    return(returnAsInput(userdata = dataset, result = Autoscore,
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

    upper <- (sum(Signature_up$Gene_Symbol %in% row.names(datasetm))/
                  nrow(Signature_up))*100
    downper <- (sum(Signature_down$Gene_Symbol %in% row.names(datasetm))/
                    nrow(Signature_down))*100
    cat(paste0("ECMSign function is using ", round(upper),
               "% of up-genes\n", "ECMSign function is using ",
               round(downper), "% of down-genes\n"))

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

    HRDSper <- (sum(HRDSdata$Gene_Symbol %in% row.names(datasetm))/
                    length(HRDSdata$Gene_Symbol))*100
    cat(paste0("HRDSSign function is using ", round(HRDSper), "% of genes\n"))

    HRDS_P <- datasetm[intersect(
        row.names(datasetm), HRDSdata[
            HRDSdata$correlation == 1, ]$Gene_Symbol), ]
    HRDS_N <- datasetm[intersect(
        row.names(datasetm), HRDSdata[
            HRDSdata$correlation == -1, ]$Gene_Symbol), ]

    HRDSscore <- unlist(lapply(X = colnames(datasetm), FUN = function(x){
                tmp <- t.test(HRDS_P[,x], HRDS_N[,x], alternative = "two.sided")
                tmp[["statistic"]]
            }))

    return(returnAsInput(userdata = dataset, result = HRDSscore,
                         SignName = "HRDS", datasetm))
}


#' Cytotoxic Immuno Signature (Cytotoxic immune cell infiltrates)
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom matrixStats colMedians
#' @importFrom stats setNames
#'
#' @import org.Hs.eg.db
#'
#' @export
CyISign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "CyISign")

    if(nametype!="SYMBOL"){
        CyIdata<- mapIds(org.Hs.eg.db, keys = CyIdata, column = nametype,
                            keytype = "SYMBOL", multiVals = "first")}
    datasetm <- getMatrix(dataset)

    ImmunoScore <- statScore(CyIdata, datasetm = datasetm, nametype = "SYMBOL",
                             typeofstat = "mean", namesignature= "CyISign")

    return(returnAsInput(userdata = dataset, result = ImmunoScore,
                         SignName = "CytoImmuno", datasetm))
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

    Lgr5per<-(sum(ISCdata$Lgr5%in%row.names(datasetm))/length(ISCdata$Lgr5))*100
    Ephper<-(sum(ISCdata$Eph%in%row.names(datasetm))/length(ISCdata$Eph))*100
    TAper<-(sum(ISCdata$TA%in%row.names(datasetm))/length(ISCdata$TA))*100
    proper<-(sum(ISCdata$pro%in%row.names(datasetm))/length(ISCdata$pro))*100

    cat(paste0("ISCSign function is using ", round(Lgr5per),
               "% of Lgr5_ISC-genes, ", round(Ephper),
               "% of EphB2_ISC-genes, ", round(TAper), "% of Late_TA-genes ",
               "and ", round(proper), "% of Proliferation-genes\n"))

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
VEGFSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "VEGFSign")

    datasetm <- getMatrix(dataset)

    VEGFscore <- log2(
        statScore(VEGFdata, datasetm = datasetm, nametype = nametype,
                  typeofstat = "mean", namesignature = "VEGFSign"))

    return(returnAsInput(userdata = dataset, result = VEGFscore,
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

    Angioscore <-statScore(
        Angiodata, datasetm = datasetm, nametype = "SYMBOL",
        typeofstat = "median", namesignature= "angioSign")

    return(returnAsInput(userdata = dataset, result = as.vector(
        scale(Angioscore)), SignName = "Angiogenesis", datasetm))
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
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
DNArepSign <- function(dataset, nametype = "SYMBOL",
                       tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "DNArepSign")

    if(nametype!="SYMBOL"){
        DNArepdata$DNAre <- mapIds(org.Hs.eg.db, keys = DNArepdata$DNAre,
                                   column = nametype, keytype = "SYMBOL",
                                   multiVals = "first")}

    datasetm <- getMatrix(dataset)

    DNAdatahigh <- DNArepdata[DNArepdata$status=="high",]
    DNAdatalow <- DNArepdata[DNArepdata$status=="low",]

    medianexp <- apply(datasetm, 2, median)

    DNArepscore <- apply(
        norm_mtx[DNAdatahigh$DNAre, ] > medianexp, 2, sum,na.rm = TRUE)
    +               apply(
        norm_mtx[DNAdatalow$DNAre, ] > medianexp, 2, sum, na.rm = TRUE)

    return(returnAsInput(userdata = dataset, result = DNArepscore,
                         SignName = "DNArepair", datasetm))
}
