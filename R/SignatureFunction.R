
#' Endothelial-Mesenchymal Transition Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the results of the Endothelial score and Mesenchymal score will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
EMTSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary", pvalues = FALSE, nperm = 100, ...) {

    firstCheck(nametype, tumorTissue, "EMTSign")

    if(nametype!="SYMBOL"){
        EMTdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = EMTdata$Gene_Symbol, column = nametype,
                                      keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    Signature_EL <- EMTdata[grep('Epithelial', EMTdata$Category),]
    Signature_ML <- EMTdata[-grep('Epithelial', EMTdata$Category),]

    cat(paste0("The function is using ", sum(Signature_EL$Gene_Symbol %in% row.names(datasetm)),
               " epithelial-like genes out of ", nrow(Signature_EL), "\nThe function is using ",
               sum(Signature_ML$Gene_Symbol %in% row.names(datasetm)), " mesenchymal-like genes out of ",
               nrow(Signature_ML),"\n"))

    gene_sets <- list(Epithelial=Signature_EL$Gene_Symbol, Mesenchymal=Signature_ML$Gene_Symbol)

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = gene_sets, method = "gsva",
                                      kcdf = "Gaussian", min.sz = 5, ssgsea.norm = FALSE, verbose = F))
    gsva_matrix <- suppressWarnings(do.call(gsva, args))

    if(pvalues){
        gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets, gsvaResult = gsva_matrix,
                                 nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(userdata = dataset, result = gsva_matrix, SignName = "", datasetm))
}


#' Pyroptosis Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Pyroptosis scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
pyroptosisSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "pyroptosisSign")

    Pyroptosisdata <- get(paste0("Pyroptosis", tumorTissue))

    if(nametype!="SYMBOL"){
        Pyroptosisdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = Pyroptosisdata$Gene_Symbol,
                                       column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    nSigGenes <- length(Pyroptosisdata$Gene_Symbol)
    cat(paste0("The function is using ", sum(Pyroptosisdata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", nSigGenes, "\n"))
    Pyroptosisdata <- Pyroptosisdata[Pyroptosisdata$Gene_Symbol %in% row.names(datasetm), ]
    Piroscore <- sapply(colnames(datasetm), function(x){
        ssgenes <- datasetm[Pyroptosisdata$Gene_Symbol, x]
        if(sum(ssgenes==0)>nSigGenes*0.5){NA}else{sum(ssgenes*Pyroptosisdata$Coefficient)}})

    return(returnAsInput(userdata = dataset, result = Piroscore, SignName = paste0("Pyroptosis",tumorTissue), datasetm))
}


#' Ferroptosis Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Ferroptosis score will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ferroptosisSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "ferroptosisSign")

    if(nametype!="SYMBOL"){
        Ferroptosisdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = Ferroptosisdata$Gene_Symbol,
                                              column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(Ferroptosisdata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", length(Ferroptosisdata$Gene_Symbol), "\n"))
    Ferroptosisdata <- Ferroptosisdata[Ferroptosisdata$Gene_Symbol %in% row.names(datasetm), ]
    ferrscore <- colSums(datasetm[Ferroptosisdata$Gene_Symbol, ]*Ferroptosisdata$Coefficient)
    return(returnAsInput(userdata = dataset, result = ferrscore, SignName = "Ferroptosis", datasetm))
}


#' Lipid Metabolism Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Lipid scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
lipidMetabolismSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary") {

    firstCheck(nametype, tumorTissue, "lipidMetabolismSign")

    if(nametype!="SYMBOL"){
        LipidMetabolismdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = LipidMetabolismdata$Gene_Symbol,
                                                  column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(LipidMetabolismdata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", length(LipidMetabolismdata$Gene_Symbol), "\n"))
    LipidMetabolismdata <- LipidMetabolismdata[LipidMetabolismdata$Gene_Symbol %in% row.names(datasetm), ]
    lipidscore <- colSums(datasetm[LipidMetabolismdata$Gene_Symbol, ] * LipidMetabolismdata$Coefficient)
    return(returnAsInput(userdata = dataset, result = lipidscore, SignName = "LipidMetabolism", datasetm))
}


#' Hypoxia Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Hypoxia scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom matrixStats colMedians
#' @importFrom stats setNames
#'
#' @import org.Hs.eg.db
#'
#' @export
hypoxiaSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue"){

    firstCheck(nametype, tumorTissue, "hypoxiaSign")

    if(nametype=="SYMBOL") { genetouse <- Hypoxiadata$Gene_Symbol
    } else if(nametype=="ENSEMBL") { genetouse <- Hypoxiadata$Gene_Ensembl
    } else (genetouse <- mapIds(org.Hs.eg.db, keys = Hypoxiadata$Gene_Symbol,
                                column = nametype, keytype = "SYMBOL", multiVals = "first"))

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(genetouse %in% rownames(datasetm)),
               " genes out of ", length(Hypoxiadata$Gene_Symbol), "\n"))
    datasetm <- datasetm[rownames(datasetm) %in% genetouse, ]

    med_counts <- colMedians(as.matrix(datasetm))

    return(returnAsInput(userdata = dataset, result = as.vector(scale(med_counts)), SignName = "Hypoxia", datasetm))
}


#' Platinum Resistance Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}} function.
#'
#' @return A SummarizedExperiment object in which the Platinum Resistance scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
platinumResistanceSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary", pvalues = FALSE, nperm = 100, ...){

    firstCheck(nametype, tumorTissue, "platinumResistanceSign")

    if(nametype!= "SYMBOL"){
        PlatinumResistancedata <- lapply(PlatinumResistancedata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype, keytype="SYMBOL", multiVals="first")))}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(PlatinumResistancedata$PlatinumResistanceUp %in% row.names(datasetm)),
              "up-genes out of", length(PlatinumResistancedata$PlatinumResistanceUp), "\nThe function is using",
              sum(PlatinumResistancedata$PlatinumResistanceDown %in% row.names(datasetm)),
              "down-genes out of", length(PlatinumResistancedata$PlatinumResistanceDown),"\n"))

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = PlatinumResistancedata, method = "gsva",
                                      kcdf = "Gaussian", min.sz = 5, ssgsea.norm = FALSE, verbose = F))
    gsva_count <- suppressWarnings(do.call(gsva, args))
    rownames(gsva_count) <- c("PlatinumResistanceUp", "PlatinumResistanceDown")

    if(pvalues){
        gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = PlatinumResistancedata,
                                 gsvaResult = gsva_matrix, nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(userdata = dataset, result = gsva_count, SignName = "", datasetm))
}


#' Prognostic high-grade serous ovarian cancer Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param age a vector of patient's age.
#' @param stage a vector of patient's tumor stage (FIGO).
#'
#' @return A SummarizedExperiment object in which the Prognostic scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
prognosticSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary", age, stage){

    firstCheck(nametype, tumorTissue, "prognosticSign")

    if(class(age)!="numeric"){stop("The age parameter must be a numeric vector")}
    if(class(stage) != "character"){stop("The stage parameter must be a character vector")}

    if(nametype!="SYMBOL"){
        names(Prognosticdata$Genes) <- mapIds(org.Hs.eg.db, keys = names(Prognosticdata$Genes),
                                        column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(names(Prognosticdata$Genes) %in% row.names(datasetm)),
              "genes out of", length(Prognosticdata$Genes), "\n"))

    intergene <- intersect(row.names(datasetm), names(Prognosticdata$Genes))
    datasetm <- datasetm[intergene,]
    gene_coeff <- colSums(datasetm*Prognosticdata$Genes[intergene])

    age_coef <- sapply(age, function(p)
        if(is.na(p)){NA
        } else {
            if(p<=53){0} else if(p>53 & p<=60){Prognosticdata$Age[1]
            } else if(p>60 & p<=67){Prognosticdata$Age[2]} else {Prognosticdata$Age[3]}})

    stage_coef <- sapply(stage, function(p)
        if(is.na(p)){Prognosticdata$Stage[2]} else if(p=="I"|p=="II"){Prognosticdata$Stage[1]} else {0})

    prog_sign <- gene_coeff+age_coef+stage_coef

    quantile_prog <- sapply(prog_sign, function(p)
        if(p<=-0.732) {"Q1"} else if(p>-0.732 & p<=-0.3126) {"Q2"
        } else if(p>-0.3126 & p<=0.0255) {"Q3"} else if(p>0.0255 & p<=0.2658) {"Q4"
        } else {"Q5"})

    return(returnAsInput(userdata = dataset, result = prog_sign, SignName = "Prognostic", datasetm))
}


#' Metabolic Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param DEdata A matrix of differentially expressed genes where rows correspond to genes,
#' the first column to Log2FoldChange and second column to its adjusted p-value,
#' or a \linkS4class{SummarizedExperiment} containing an assay in which rows correspond to genes, the first column to Log2FoldChange and second column to its adjusted p-value.
#' @param nametype gene name ID of your DEdata (row names).
#' @param nsamples number of samples in the DEdata.
#'
#' @return A data frame with a Metabolic score for each pathway and the respective p-values.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
metabolicSign <- function(DEdata, nametype = "SYMBOL", tumorTissue = "pan-tissue", nsamples){

    firstCheck(nametype, tumorTissue, "metabolicSign")

    if(class(nsamples)!="numeric"){stop("The nsample parameter must be a numeric vector")}

    gene_score <- abs(DEdata[,1] -log(DEdata[,2]))
    names(gene_score) <- row.names(DEdata)

    if(nametype!="SYMBOL"){
        Metabolismdata <- lapply(Metabolismdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype, keytype="SYMBOL", multiVals="first")))}

    gene_pathway <- lapply(Metabolismdata, intersect, row.names(DEdata))
    path_score <- sapply(gene_pathway, function(x) sum(gene_score[x]))/sqrt(nsamples)
    pvals <- c()
    for(i in 1:length(path_score)){
        z <- c()
        for(j in 1:10000){
            bootscore <- sample(gene_score, size=lengths(gene_pathway)[i], replace = T)
            z[j] <- sum(bootscore)/sqrt(nsamples)}
        pvals[i] <- sum(z>=path_score[i])/10000
    }
    return(cbind(MetabolicScore=path_score, Pvalue=pvals))
}


#' Immunogenic Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Immunogenic scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
immunoScoreSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary", author = "Hao"){

    firstCheck(nametype, tumorTissue, "immunoScoreSign")

    if(nametype!="SYMBOL"){
        ImmunoScoredata$genes <- mapIds(org.Hs.eg.db, keys = ImmunoScoredata$genes, column = nametype,
                                    keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    g <- intersect(row.names(datasetm), ImmunoScoredata$genes)

    cat(paste("The function is using", length(g), "genes out of", length(ImmunoScoredata$genes), "\n"))

    subdataset <- datasetm[g,]
    ImmunoScoredata <- ImmunoScoredata[ImmunoScoredata$genes %in% g, ]

    SE <- (ImmunoScoredata$HR - ImmunoScoredata$`95CI_L`)/1.96
    k <- (1 - ImmunoScoredata$HR)/SE

    ImmunoScores <- unlist(lapply(seq_len(ncol(subdataset)), function(p) sum(k*subdataset[,p], na.rm = T)))

    return(returnAsInput(userdata = dataset, result = ImmunoScores, SignName = "ImmunoScore", datasetm))
}


#' ConsensusOV Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param method the subtyping method to use. Default is "consensusOV".
#' @param ... optional parameters to be passed to \code{\link[consensusOV]{get.subtypes}}.
#'
#' @return A SummarizedExperiment object in which the COnsensusOV scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom consensusOV get.subtypes
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
consensusOVSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary", method = "consensusOV", ...){

    firstCheck(nametype, tumorTissue, "consensusOVSign")

    datasetm <- getMatrix(dataset)

    if(nametype!="ENTREZID"){
        genename <- mapIds(org.Hs.eg.db, keys = row.names(datasetm), column = "ENTREZID",
                           keytype = nametype, multiVals = "first")
        datasetm <- datasetm[!is.na(genename),]
        genename <- genename[!is.na(genename)]
        datasetm <- datasetm[!duplicated(genename),]
        genename <- genename[!duplicated(genename)]
    } else {genename <- row.names(datasetm)}

    consensus_subtypes <- get.subtypes(expression.dataset=datasetm, entrez.ids=genename, method=method, ...)

    return(returnAsInput(userdata = dataset, result = t(consensus_subtypes$rf.probs), SignName = "", datasetm))
}


#' ImmunoPhenoScore Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the IPS, MHC, CP, EC and SC scores will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
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
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype, keytype="SYMBOL", multiVals="first"))))}

    datasetm <- getMatrix(dataset)
    sample_names <- colnames(datasetm)

    cat(paste("The function is using", sum(IPSdata$GENE %in% row.names(datasetm)),
              "genes out of", nrow(IPSdata), "\n"))

    ind <- which(is.na(match(IPSdata$GENE, row.names(datasetm))))
    MISSING_GENES <- IPSdata$GENE[ind]
    if (length(MISSING_GENES)>0) {cat("Differently named or missing genes: ", MISSING_GENES, "\n")}

    IPS <- NULL; MHC <- NULL; CP <- NULL; EC <- NULL; SC <- NULL; AZ <- NULL
    for (i in 1:length(sample_names)) {
        GE <- datasetm[,i]
        Z1 <- (datasetm[match(IPSdata$GENE, row.names(datasetm)),i]-mean(GE))/sd(GE)
        WEIGHT <- NULL; MIG <- NULL; k <- 1
        for (gen in unique(IPSdata$NAME)) {
            MIG[k] <- mean(Z1[which(IPSdata$NAME==gen)], na.rm=TRUE)
            WEIGHT[k] <- mean(IPSdata$WEIGHT[which(IPSdata$NAME==gen)])
            k<-k+1}
        WG <- MIG*WEIGHT
        MHC[i]<-mean(WG[1:10], na.rm = T)
        CP[i]<-mean(WG[11:20], na.rm = T)
        EC[i]<-mean(WG[21:24], na.rm = T)
        SC[i]<-mean(WG[25:26], na.rm = T)
        AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
        IPS[i]<-ipsmap(AZ[i])}

    ipsres <- data.frame(IPS, MHC, CP, EC, SC)
    row.names(ipsres) <- sample_names
    return(returnAsInput(userdata = dataset, result = t(ipsres), SignName = "", datasetm))
}


#' Core Matrisome Gene signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the median gene expression based on the core matrisome signature will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom matrixStats colMedians
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
matrisomeSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue") {

    firstCheck(nametype, tumorTissue, "matrisomeSign")

    if(nametype!="SYMBOL"){
        Matrisomedata <- mapIds(org.Hs.eg.db, keys=Matrisomedata, column=nametype,
                                keytype="SYMBOL", multiVals="first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(Matrisomedata %in% row.names(datasetm)), "matrisome's genes out of 9\n"))

    median_cm <- colMedians(datasetm[row.names(datasetm) %in% Matrisomedata, ])

    return(returnAsInput(userdata = dataset, result = median_cm, SignName = "Matrisome", datasetm))
}

#' Mitotic Index
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the means gene expression based on the mitotix index will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
mitoticIndexSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "pan-tissue") {

    firstCheck(nametype, tumorTissue, "mitoticIndexSign")

    if(nametype!="SYMBOL"){
        MitoticIndexdata <- mapIds(org.Hs.eg.db, keys=MitoticIndexdata, column=nametype,
                                   keytype="SYMBOL", multiVals="first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(MitoticIndexdata %in% row.names(datasetm)),
              "mititotic index genes out of 9\n"))

    MI_means <- colMeans(datasetm[row.names(datasetm) %in% MitoticIndexdata, ])

    return(returnAsInput(userdata = dataset, result = MI_means, SignName = "MitoticIndex", datasetm))
}

#' Local Immune Cytolytic Activity (CYT) Signature
#'
#' This signature is computed accordingly to the reference paper,
#' to have more details explore the function \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Expression values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' Alternatively an object of type \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \linkS4class{SpatialExperiment} or \linkS4class{Seurat}
#' containing an assay where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param tumorTissue tissue of the tumor.
#'
#' @return A SummarizedExperiment object in which the means gene expression based on the mitotix index will be added
#' in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom labstatR meang
#' @import org.Hs.eg.db
#'
#' @export
CYTSign <- function(dataset, nametype = "SYMBOL", tumorTissue = "ovary"){

    firstCheck(nametype, tumorTissue, "CYTSign")

    if(nametype!="SYMBOL"){
        CYTdata <- mapIds(org.Hs.eg.db, keys = CYTdata, column = nametype,
                           keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(CYTdata %in% row.names(datasetm)),
              "genes out of", length(CYTdata), "\n"))

    CYTdata <- CYTdata[CYTdata %in% row.names(datasetm)]
    CYTScore <- apply(datasetm[CYTdata,], 2, meang)

    return(returnAsInput(userdata = dataset, result = CYTScore, SignName = "CYT", datasetm))
}
