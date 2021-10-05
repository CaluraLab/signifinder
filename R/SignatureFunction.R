
#' Endothelial-Mesenchymal Transition Signature
#'
#' Given a dataset, EMTSign returns the Endothelial score and the Mesenchymal score for
#' each sample, based on the work of QH Miow at al. (2015).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return A SummarizedExperiment object in which the results of the Endothelial score and Mesenchymal score will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Miow, Q., Tan, T., Ye, J. et al. Epithelial–mesenchymal status renders differential responses to cisplatin in ovarian cancer.
#' Oncogene 34, 1899–1907 (2015). \url{https://doi.org/10.1038/onc.2014.136}
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


#' Pyroptosis Signature
#'
#' Given a dataset, pyroptosisSign returns the pyroptosis score for each sample, based on Ying Ye et al. (2021).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Pyroptosis scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Ye Y, Dai Q. & Qi H. A novel defined pyroptosis-related gene signature for predicting the prognosis of ovarian cancer.
#' Cell Death Discov. 7, 71 (2021). \url{https://doi.org/10.1038/s41420-021-00451-x}
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


#' Ferroptosis Signature
#'
#' Given a dataset, ferroptosisSign returns the Ferroptosis score for each sample Ying Ye et al. (2021).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Ferroptosis score will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Ye Y, Dai Q, Li S, He J and Qi H (2021) A Novel Defined Risk Signature of the Ferroptosis-Related Genes for Predicting the Prognosis of Ovarian Cancer.
#' Front. Mol. Biosci. 8:645845. \url{https://doi.org/10.3389/fmolb.2021.645845}
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


#' Lipid Metabolism Signature
#'
#' Given a dataset, lipidMetabolismSign returns the Lipid score for each sample Mingjun Zheng et al. (2020).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Lipid scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Zheng M, Mullikin H, Hester A, Czogalla B, Heidegger H, Vilsmaier T, Vattai A, Chelariu-Raicu A, Jeschke U, Trillsch F, Mahner S, Kaltofen T. Development and Validation of a Novel 11-Gene Prognostic Model for Serous Ovarian Carcinomas Based on Lipid Metabolism Expression Profile.
#' International Journal of Molecular Sciences. 2020; 21(23):9169. \url{https://doi.org/10.3390/ijms21239169}
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
#' Given a dataset, Hypoxia returns the hypoxia score for each sample as in Francesca M. Buffa et al. 2010.
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Hypoxia scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Buffa F, Harris A, West C. et al. Large meta-analysis of multiple cancers reveals a common, compact and highly prognostic hypoxia metagene.
#' Br J Cancer 102, 428–435 (2010). \url{https://doi.org/10.1038/sj.bjc.6605450}
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
#' Given a dataset, it returns the gsva scores for each sample from International Cancer Genome Consortium (ICGC).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param pvalues logical. It allows to compute p-values by permutations.
#' @param nperm number of permutations.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return A SummarizedExperiment object in which the Platinum Resistance scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
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


#' Prognostic high-grade serous ovarian cancer Signature
#'
#' Given a dataset, based on a 101-gene expression signature by J. Millstein et. al (2020), prognosticSign returns the Quantile assignation for each sample
#' allowing to improve risk stratification in clinical trials by identifying patients who are least likely to achieve 5-year survival.
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param age a vector of patient's age.
#' @param stage a vector of patient's tumor stage (FIGO).
#'
#' @return A SummarizedExperiment object in which the Prognostic scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Millstein J, Budden T, Goode EL, Anglesio MS, Talhouk A, Intermaggio MP, et al. Prognostic gene expression signature for high-grade serous ovarian cancer.
#' Ann Oncol. 2020;31(9):1240–50. \url{https://doi.org/10.1016/j.annonc.2020.05.019}
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
#' Given a list of DEG, metabolicSign returns a matrix with pathways score and a correspondent p-value calculated with Bootstrapping.
#' The signature is based on the work of Rosario et. al (2018).
#'
#' @param DEdata A matrix of differentially expressed genes where rows correspond to genes,
#' the first column to Log2FoldChange and second column to its adjusted p-value,
#' or a SummarizedExperiment which contains an assay represented by a matrix-like object of differentially expressed genes
#' where rows correspond to genes, the first column to Log2FoldChange and second column to its adjusted p-value.
#' @param nametype gene name ID of your DEdata (row names).
#' @param nsamples number of samples in the DEdata.
#'
#' @return A SummarizedExperiment object in which the Metabolic scores and their respective p-values will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Rosario SR, Long MD, Affronti H.C. et al. Pan-cancer analysis of transcriptional metabolic dysregulation using The Cancer Genome Atlas.
#' Nat Commun 9, 5330 (2018). \url{https://doi.org/10.1038/s41467-018-07232-8}
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
#' Given a dataset, immunoScoreSign returns the ImmunoScore for each sample. This signature is
#' based on Dapeng Hao et. al (2018).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the Immunogenic scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Hao D, Liu J, Chen M, Li J, Wang L, Li X, Zhao Q, Di LJ. Immunogenomic Analyses of Advanced Serous Ovarian Cancer Reveal Immune Score is a Strong Prognostic Factor and an Indicator of Chemosensitivity.
#' Clin Cancer Res. 2018 Aug 1;24(15):3560-3571. \url{https://doi.org/10.1158/1078-0432.CCR-17-3862}
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
#' Given a dataset, consensusOVSign returns ovarian cancer subtypes. This signature is based on Chen et. al (2018).
#'
#' @param dataset Expression matrix. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#' @param method the subtyping method to use. Default is "consensusOV".
#' @param ... optional parameters to be passed to the low level function.
#'
#' @return A SummarizedExperiment object in which the COnsensusOV scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Gregory M. Chen, Lavanya Kannan, Ludwig Geistlinger, Victor Kofia et al. Consensus on Molecular Subtypes of High-Grade Serous Ovarian Carcinoma
#' Clin Cancer Res October 15 2018 (24) (20) 5037-5047. \url{https://doi.org/10.1158/1078-0432.CCR-18-0784}
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
#' Given a dataset, IPSSign returns fro each sample, the IPSs (ImmunoPhenoScores) and scores for the four groups of immune biomarkers: MHC (Antigen processing molecules),
#' CP (Checkpoints immunomodulators), EC (Immune effector cells) and SC (Suppressor cells) identified by P Charoentong et. al (2017).
#'
#' @param dataset Expression matrix of TPM values. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the IPS, MHC, CP, EC and SC scores will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Charoentong P, Finotello F, Angelova M, Mayer C, Efremova M, Rieder D, Hackl H, Trajanoski Z. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade.
#' Cell Rep. 2017 Jan 3;18(1):248-262. \url{https://doi.org/10.1016/j.celrep.2016.12.019}
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
#' Given a dataset, matrisomeSign returns the median genes expression based on Yuzhalin et all. (2018).
#'
#' @param dataset Expression matrix in which row names must be Official Symbol. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the median gene expression based on the core matrisome signature will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Yuzhalin A, Urbonas T, Silva M. et al. A core matrisome gene signature predicts cancer outcome.
#' British Journal of Cancer volume 118, 435–440 (2018). \url{https://doi.org/10.1038/bjc.2017.458}
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
#' Given a dataset, mitoticIndexSign returns the means genes expression based on Yang et all. (2016).
#'
#' @param dataset Expression matrix in which row names must be Official Symbol. A data frame or a matrix where rows correspond to genes and columns correspond to samples.
#' A SummarizedExperiment which contains an assay represented by a matrix-like object of numeric. The rows typically represent
#' genomic ranges of interest and the columns represent samples.
#' @param nametype gene name ID of your dataset (row names).
#'
#' @return A SummarizedExperiment object in which the means gene expression based on the mitotix index will be added
#' in the `colData` section which contains sample meta-data describing the samples.
#'
#' @references Yang, Z., Wong, A., Kuh, D. et al. Correlation of an epigenetic mitotic clock with cancer risk.
#' Genome Biol 17, 205 (2016). \url{https://doi.org/10.1186/s13059-016-1064-3}
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
