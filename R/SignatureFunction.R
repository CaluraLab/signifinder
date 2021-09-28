
#' Endothelial-Mesenchymal Transition Signature
#'
#' Given a dataset, it returns the Endothelial score and the Mesenchymal score for
#' each sample, based on QH Miow at al. (2015).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param pvalues whether to compute pvalues by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
EMTSign <- function(dataset, nametype = "SYMBOL", pvalues = FALSE, nperm = 100, ...) {

    firstCheck(nametype)

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


#' Piroptosis Signature
#'
#' Given a dataset, it returns the piroptosis score for each sample, based on Mingjun Zheng et al. (2020).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
piroptosisSign <- function(dataset, nametype = "SYMBOL"){

    firstCheck(nametype)

    if(nametype!="SYMBOL"){
        Piroptosisdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = Piroptosisdata$Gene_Symbol,
                                       column = nametype, keytype = "SYMBOL", multiVals = "first")}

    datasetm <- getMatrix(dataset)

    nSigGenes <- length(Piroptosisdata$Gene_Symbol)
    cat(paste0("The function is using ", sum(Piroptosisdata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", nSigGenes, "\n"))
    Piroptosisdata <- Piroptosisdata[Piroptosisdata$Gene_Symbol %in% row.names(datasetm), ]
    Piroscore <- sapply(colnames(datasetm), function(x){
        ssgenes <- datasetm[Piroptosisdata$Gene_Symbol, x]
        if(sum(ssgenes==0)>nSigGenes*0.5){NA}else{sum(ssgenes*Piroptosisdata$Coefficient)}})
    # Piroscore <- colSums(datasetm[Piroptosisdata$Gene_Symbol, ]*Piroptosisdata$Coefficient)
    return(returnAsInput(userdata = dataset, result = Piroscore, SignName = "Piroptosis", datasetm))
}


#' FerroptosisSignature
#'
#' Given a dataset, it returns the Ferroptosis score for each sample Ying Ye et al. (2021).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
ferroptosisSign <- function(dataset, nametype = "SYMBOL"){

    firstCheck(nametype)

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


#' LipidSignature
#'
#' Given a dataset, it returns the Lipid score for each sample Mingjun Zheng et al. (2020).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
lipidMetabolismSign <- function(dataset, nametype = "SYMBOL") {

    firstCheck(nametype)

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
#' Given a dataset, it returns the hypoxia score for each sample as in Buffa et al. 2010.
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom matrixStats colMedians
#' @importFrom stats setNames
#'
#' @import org.Hs.eg.db
#'
#' @export
hypoxiaSign <- function(dataset, nametype = "SYMBOL"){

    firstCheck(nametype)

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
#' Given a dataset, it returns the gsva score for each sample from International Cancer Genome Consortium (ICGC).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param pvalues whether to compute pvalues by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
platinumResistanceSign <- function(dataset, nametype = "SYMBOL", pvalues = FALSE, nperm = 100, ...){

    firstCheck(nametype)

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


#' Prognostic Signature
#'
#' Given a dataset, it returns the Quantile assignation for each sample from J. Millstein et. al (2020).
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param age vector of patient's age.
#' @param stage vector of patient's tumor stage (FIGO).
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
prognosticSign <- function(dataset, nametype = "SYMBOL", age, stage){

    firstCheck(nametype)

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
        if(p<=53){0} else if(p>53 & p<=60){Prognosticdata$Age[1]
        } else if(p>60 & p<=67){Prognosticdata$Age[2]} else {Prognosticdata$Age[3]})

    stage_coef <- sapply(stage, function(p)
        if(p=="NA"){Prognosticdata$Stage[2]} else if(p=="I"|p=="II"){Prognosticdata$Stage[1]} else {0})

    prog_sign <- gene_coeff+age_coef+stage_coef

    quantile_prog <- sapply(prog_sign, function(p)
        if(p<=-0.732) {"Q1"} else if(p>-0.732 & p<=-0.3126) {"Q2"
        } else if(p>-0.3126 & p<=0.0255) {"Q3"} else if(p>0.0255 & p<=0.2658) {"Q4"
        } else {"Q5"})

    return(returnAsInput(userdata = dataset, result = prog_sign, SignName = "Prognostic", datasetm))
}


#' Metabolic Signature
#'
#' Given a list of DEG, it returns a matrix with pathways score and a correspondent pvalue
#' calculated with Bootstrapping.
#' This signature is based on Rosario et. al (2018).
#'
#' @param DEdata matrix of differentially expressed genes where rows correspond to genes,
#' first column to Log2FoldChange and second column to its adjusted pvalue.
#' @param nametype gene name ID of your DEdata row names.
#' @param nsamples number of samples in the DEdata.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
metabolicSign <- function(DEdata, nametype = "SYMBOL", nsamples){

    firstCheck(nametype)

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
#' Given a dataset, it returns the ImmunoScore for each sample. This signature is
#' based on Dapeng Hao et. al (2018).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
immunoScoreSign <- function(dataset, nametype = "SYMBOL"){

    firstCheck(nametype)

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
#' Given a dataset, it returns ovarian cancer subtypes. This signature is based on Chen et. al (2018).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param method the subtyping method to use. Default is "consensusOV".
#' @param ... optional parameters to be passed to the low level function.
#'
#' @return NULL
#'
#' @importFrom consensusOV get.subtypes
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
consensusOVSign <- function(dataset, nametype = "SYMBOL", method = "consensusOV", ...){

    firstCheck(nametype)

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
#' Given a dataset, it returns IPS for each sample. This signature is based on .. et. al (..).
#'
#' @param dataset TPM expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices pdf dev.off
#' @importFrom stats sd
#'
#' @export
IPSSign <- function(dataset, nametype = "SYMBOL"){

    firstCheck(nametype)

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
#' Given a dataset, it returns the median genes expression based on Yuzhalin et all. (2018).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' Row names must be Official Symbol.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom matrixStats colMedians
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
matrisomeSign <- function(dataset, nametype = "SYMBOL") {

    firstCheck(nametype)

    if(nametype!="SYMBOL"){
        Matrisomedata <- mapIds(org.Hs.eg.db, keys=Matrisomedata, column=nametype,
                                keytype="SYMBOL", multiVals="first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(Matrisomedata %in% row.names(dataset)), "matrisome's genes out of 9\n"))

    median_cm <- colMedians(datasetm[row.names(datasetm) %in% Matrisomedata, ])

    return(returnAsInput(userdata = dataset, result = median_cm, SignName = "Matrisome", datasetm))
}

#' Mitotic Index
#'
#' Given a dataset, it returns the means genes expression based on Yang et all. (2016).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples.
#' Row names must be Official Symbol.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
mitoticIndexSign <- function(dataset, nametype = "SYMBOL") {

    firstCheck(nametype)

    if(nametype!="SYMBOL"){
        MitoticIndexdata <- mapIds(org.Hs.eg.db, keys=MitoticIndexdata, column=nametype,
                                   keytype="SYMBOL", multiVals="first")}

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(MitoticIndexdata %in% row.names(datasetm)),
              "mititotic index genes out of 9\n"))

    MI_means <- colMeans(datasetm[row.names(datasetm) %in% MitoticIndexdata, ])

    return(returnAsInput(userdata = dataset, result = MI_means, SignName = "MitoticIndex", datasetm))
}
