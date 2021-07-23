
matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)}

getMatrix <- function(userdata){
    if(!is.matrix(userdata)) {
        if(class(userdata)=="Seurat"){
            if(length(userdata@assays)==1){userdata <- as.matrix(userdata@assays$RNA@data)
            } else {userdata <- as.matrix(userdata@assays$SCT@data)}
        } else if(class(userdata)%in%c("SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
            userdata <- as.matrix(SummarizedExperiment::assay(userdata))
        } else if(class(userdata)=="data.frame"){userdata <- as.matrix(userdata)
        } else {stop("Dataset is not supported")}}
    return(userdata)}

returnAsInput <- function(userdata, result, SignName){
    if(!is.matrix(userdata) & !is.data.frame(userdata)) {
        if(class(userdata)=="Seurat"){
            if(is.vector(result)){names <- c(colnames(userdata@meta.data), SignName)
                userdata@meta.data <- cbind(userdata@meta.data, name=result)
                colnames(userdata@meta.data) <- names
            } else {userdata@meta.data <- cbind(userdata@meta.data, t(result))}
        } else if(class(userdata)%in%c("SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
            names <- c(colnames(userdata@colData), SignName)
            if(is.vector(result)){userdata@colData <- cbind(userdata@colData, name=result)
            colnames(userdata@colData) <- names
            } else {userdata@colData <- cbind(userdata@colData, t(result))}}
        return(userdata)
    } else {return(result)}}

ipsmap <- function(x){
    if (x<=0) {ips <- 0} else if (x>=3) {ips <- 10} else {ips <- round(x*10/3, digits=0)}
    return(ips)}

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)

mapcolors <- function(x){
    za <- NULL
    if (x>=3) {za = 1000} else if (x<=-3) {za = 1} else {za = round(166.5*x+500.5, digits = 0)}
    return(my_palette[za])}

mapbw <- function(x){
    za <- NULL
    if (x>=2) {za = 1000} else if (x<=-2) {za = 1} else {za = round(249.75*x+500.5, digits = 0)}
    return(my_palette2[za])}

GSVAPvalues <- function(expr, gset.idx.list, gsvaResult, nperm, args){
    datasetGenes <- rownames(expr)
    filteredGeneSets <- lapply(gset.idx.list, y = datasetGenes, intersect)
    permutedResults <- lapply(seq_len(nperm), function(x){
        permlist <- lapply(seq_len(length(gset.idx.list)), function(i)
            sample(datasetGenes, size = lengths(filteredGeneSets)[i], replace = F))
        gsva_matrix <- suppressWarnings(do.call(gsva, args))
        data.frame(t(gsva_matrix))})
    permutedResByGeneSet <- split.default(x = Reduce(cbind, permutedResults), seq_len(length(gset.idx.list)))
    permutedResByGeneSet <- lapply(permutedResByGeneSet, function(x)data.frame(t(x)))
    finalRes <- do.call(rbind, lapply(seq_len(length(gset.idx.list)), function(i){
        gspvalues <- sapply(1:ncol(expr), function(j){
            (min(c(sum(permutedResByGeneSet[[i]][,j]<=gsvaResult[i,j]),
                   sum(permutedResByGeneSet[[i]][,j]>=gsvaResult[i,j])))+1)/(nperm+1)})
        gspvalues}))
    colnames(finalRes) <- colnames(expr)
    rownames(finalRes) <- paste(names(gset.idx.list), "pval", sep = "_")
    return(finalRes)}


#' 7 spots
#'
#' Given a 10X Visium dataset, it reassigns to each spot the aggregation of it with the nearest.
#'
#' @param dataset Seurat object of a 10X Visium dataset
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
GetAggregatedSpot <- function(dataset){
    spotcoords <- data.frame(row = as.vector(sapply(seq(0,76,2), function(i) rep(c(i,i+1),64) )),
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
    ## now how to add counts to the seurat object??
    return(dataset)
}


#' Endothelial-Mesenchymal Transition Signature
#'
#' Given a dataset, it returns the Endothelial score and the Mesenchymal score for each sample, based on QH Miow at al. (2015).
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
EMTSign <- function(dataset, nametype, pvalues = FALSE, nperm = 100, ...) {

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        EMTdata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys=EMTdata$Gene_Symbol, column=nametype, keytype="SYMBOL", multiVals="first")
    }

    datasetm <- getMatrix(dataset)

    Signature_EL <- EMTdata[grep('Epithelial-like', EMTdata$Category),]
    Signature_ML <- EMTdata[-grep('Epithelial-like', EMTdata$Category),]

    cat(paste0("The function is using ", sum(Signature_EL$Gene_Symbol %in% row.names(datasetm)),
               " epithelial-like genes out of ", nrow(Signature_EL), "\nThe function is using ",
               sum(Signature_ML$Gene_Symbol %in% row.names(datasetm)), " mesenchymal-like genes out of ",
               nrow(Signature_ML),"\n"))

    gene_sets <- list(Epithelial=Signature_EL$Gene_Symbol, Mesenchymal=Signature_ML$Gene_Symbol)

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = gene_sets, method = "gsva",
                                      kcdf = "Gaussian", min.sz = 5, ssgsea.norm = FALSE))
    gsva_matrix <- do.call(gsva, args)

    if(pvalues){
        gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets, gsvaResult = gsva_matrix,
                                 nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(userdata = dataset, result = gsva_matrix, SignName = ""))
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
PiroSign <- function(dataset, nametype){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        Pirodata$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = Pirodata$Gene_Symbol,
                                       column = nametype, keytype = "SYMBOL", multiVals = "first")
    }

    datasetm <- getMatrix(dataset)

    nSigGenes <- length(Pirodata$Gene_Symbol)
    cat(paste0("The function is using ", sum(Pirodata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", nSigGenes, "\n"))
    Pirodata <- Pirodata[Pirodata$Gene_Symbol %in% row.names(datasetm), ]
    Piroscore <- sapply(colnames(datasetm), function(x){
        ssgenes <- datasetm[Pirodata$Gene_Symbol, x]
        if(sum(ssgenes==0)>nSigGenes*0.5){NA}else{sum(ssgenes*Pirodata$Coefficient)}})
    # Piroscore <- colSums(datasetm[Pirodata$Gene_Symbol, ]*Pirodata$Coefficient)
    return(returnAsInput(userdata = dataset, result = Piroscore, SignName = "PiroptosisScore"))
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
#' @param nametype gene name ID of rownames of dataset.
#' @import org.Hs.eg.db
#'
#' @export
FerrSign <- function(dataset, nametype){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        Ferrdata$Gene_Symbol <- mapIds(org.Hs.eg.db,keys= Ferrdata$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(Ferrdata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", length(Ferrdata$Gene_Symbol), "\n"))
    Ferrdata <- Ferrdata[Ferrdata$Gene_Symbol %in% row.names(datasetm), ]
    ferrscore <- colSums(datasetm[Ferrdata$Gene_Symbol, ]*Ferrdata$Coefficient)
    return(returnAsInput(userdata = dataset, result = ferrscore, SignName = "FerroptosisScore"))
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
LipidMetSign <- function(dataset, nametype) {

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        Lipidata$Gene_Symbol <- mapIds(org.Hs.eg.db,keys= Lipidata$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(Lipidata$Gene_Symbol %in% row.names(datasetm)),
               " genes out of ", length(Lipidata$Gene_Symbol), "\n"))
    Lipidata <- Lipidata[Lipidata$Gene_Symbol %in% row.names(datasetm), ]
    lipidscore <- colSums(datasetm[Lipidata$Gene_Symbol, ] * Lipidata$Coefficient)
    return(returnAsInput(userdata = dataset, result = lipidscore, SignName = "LipidScore"))
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
HypoSign <- function(dataset, nametype){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype=="SYMBOL") { genetouse <- Hypodata$Gene_Symbol
    } else if(nametype=="ENSEMBL") { genetouse <- Hypodata$Gene_Ensembl
    } else (genetouse <- mapIds(org.Hs.eg.db,keys= Hypodata$Gene_Symbol,
                                column= nametype, keytype="SYMBOL", multiVals="first"))

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(genetouse %in% rownames(datasetm)),
               " genes out of ", length(Hypodata$Gene_Symbol), "\n"))
    datasetm <- datasetm[rownames(datasetm) %in% genetouse, ]

    med_counts <- setNames(colMedians(as.matrix(datasetm)), colnames(datasetm))
    scores <- data.frame(E = med_counts, HS = scale(med_counts))

    return(returnAsInput(userdata = dataset, result = t(scores), SignName = ""))
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
PlatResSign <- function(dataset, nametype, pvalues = FALSE, nperm = 100, ...){

    if (!nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS")){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!= "SYMBOL"){
        Platdata <- lapply(Platdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys = x, column = nametype, keytype = "SYMBOL", multiVals = "first")))
    }

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(Platdata$up %in% row.names(datasetm)),
              "up-genes out of", length(Platdata$up), "\nThe function is using",
              sum(Platdata$down %in% row.names(datasetm)), "down-genes out of", length(Platdata$down),"\n"))

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = Platdata, method = "gsva",
                                      kcdf = "Gaussian", min.sz = 5, ssgsea.norm = FALSE))
    gsva_count <- do.call(gsva, args)

    if(pvalues){
        gsva_pval <- GSVAPvalues(expr = datasetm, gset.idx.list = gene_sets, gsvaResult = gsva_matrix,
                                 nperm = nperm, args = args)
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)}

    return(returnAsInput(userdata = dataset, result = gsva_count, SignName = ""))
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
PrognosticSign <- function(dataset, nametype, age, stage){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }
    if(class(age)!="numeric"){
        stop("The age parameter must be a numeric vector")
    }
    if(class(stage) != "character"){
        stop("The stage parameter must be a character vector")
    }

    if(nametype!="SYMBOL"){
        names(Progdata$Genes) <- mapIds(org.Hs.eg.db, keys = names(Progdata$Genes),
                                        column = nametype, keytype = "SYMBOL", multiVals = "first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(names(Progdata$Genes) %in% row.names(datasetm)),
              "genes out of", length(Progdata$Genes), "\n"))

    intergene <- intersect(row.names(datasetm), names(Progdata$Genes))
    datasetm <- datasetm[intergene,]
    gene_coeff <- colSums(datasetm*Progdata$Genes[intergene])

    age_coef <- sapply(age, function(p)
        if(p<=53){0} else if(p>53 & p<=60){Progdata$Age[1]
        } else if(p>60 & p<=67){Progdata$Age[2]} else {Progdata$Age[3]})

    stage_coef <- sapply(stage, function(p)
        if(p=="NA"){Progdata$Stage[2]} else if(p=="I"|p=="II"){Progdata$Stage[1]} else {0})

    prog_sign <- gene_coeff+age_coef+stage_coef

    quantile_prog <- sapply(prog_sign, function(p)
        if(p<=-0.732) {"Q1"} else if(p>-0.732 & p<=-0.3126) {"Q2"
        } else if(p>-0.3126 & p<=0.0255) {"Q3"} else if(p>0.0255 & p<=0.2658) {"Q4"
        } else {"Q5"})

    return(returnAsInput(userdata = dataset, result = prog_sign, SignName = "PrognosticScore"))
}


#' Metabolic Signature
#'
#' Given a list of DEG, it returns a matrix with pathways score and a correspondent pvalue calculated with Bootstrapping.
#' This signature is based on Rosario et. al (2018).
#'
#' @param DEdata matrix of differentially expressed genes where rows correspond to genes, first column to Log2FoldChange and second column to its adjusted pvalue.
#' @param nametype gene name ID of your DEdata row names.
#' @param nsamples number of samples in the DEdata.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
MetabolicSign <- function(DEdata, nametype, nsamples){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }
    if(class(nsamples)!="numeric"){
        stop("The nsample parameter must be a numeric vector")
    }

    gene_score <- abs(DEdata[,1] -log(DEdata[,2]))
    names(gene_score) <- row.names(DEdata)

    if(nametype!="SYMBOL"){
        Metadata <- lapply(Metadata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }

    gene_pathway <- lapply(Metadata, intersect, row.names(DEdata))
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
#' Given a dataset, it returns the ImmunoScore for each sample. This signature is based on Dapeng Hao et. al (2018).
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
ImmunoSign <- function(dataset, nametype){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        Immudata$genes <- mapIds(org.Hs.eg.db, keys = Immudata$genes, column = nametype,
                                    keytype = "SYMBOL", multiVals = "first")
    }

    datasetm <- getMatrix(dataset)

    g <- intersect(row.names(datasetm), Immudata$genes)

    subdataset <- datasetm[g,]
    Immudata <- Immudata[Immudata$genes %in% g, ]

    SE <- (Immudata$HR - Immudata$`95CI_L`)/1.96
    k <- (1 - Immudata$HR)/SE

    ImmunoScores <- unlist(lapply(seq_len(ncol(subdataset)), function(p) sum(k*subdataset[,p], na.rm = T)))

    return(returnAsInput(userdata = dataset, result = ImmunoScores, SignName = "ImmunoScore"))
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
ConsensusOVSign <- function(dataset, nametype, method = "consensusOV", ...){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    datasetm <- getMatrix(dataset)

    if(nametype!="ENTREZID"){
        genename <- mapIds(org.Hs.eg.db, keys = row.names(datasetm), column = "ENTREZID",
                           keytype = nametype, multiVals = "first")
        datasetm <- datasetm[!is.na(genename),]
        genename <- genename[!is.na(genename)]
        datasetm <- datasetm[!duplicated(genename),]
        genename <- genename[!duplicated(genename)]
    } else {genename <- row.names(datasetm)}

    consensus_subtypes <- get.subtypes(expression.dataset = datasetm, entrez.ids = genename, method = method, ...)

    return(returnAsInput(userdata = dataset, result = consensus_subtypes, SignName = ""))
}


#' ImmunoPhenoScore Signature
#'
#' Given a dataset, it returns IPS for each sample. This signature is based on .. et. al (..).
#'
#' @param dataset TPM expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param makeplot whether to make or not IPS plots. These will be one for each sample. Default is "no".
#' @param plotpath directory path to save plots.
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
IPSSign <- function(dataset, nametype, makeplot = "no", plotpath=NULL){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(makeplot=="yes"){if(is.null(plotpath)){
        message("Plots will be saved in the working directory because no path has been specified.")}
    }

    if(nametype!="SYMBOL"){
        IPSGdata[,c(1,2)] <- data.frame(lapply(IPSGdata[,c(1,2)], function(x)
            suppressMessages(mapIds(org.Hs.eg.db, keys=x, column=nametype, keytype="SYMBOL", multiVals="first"))))
    }

    datasetm <- getMatrix(dataset)
    sample_names <- colnames(datasetm)

    cat(paste("The function is using", sum(IPSGdata$GENE %in% row.names(datasetm)),
              "genes out of ", nrow(IPSGdata), "\n"))

    ind <- which(is.na(match(IPSGdata$GENE, row.names(datasetm))))
    MISSING_GENES <- IPSGdata$GENE[ind]
    if (length(MISSING_GENES)>0) {cat("Differently named or missing genes: ", MISSING_GENES, "\n")}

    IPS <- NULL; MHC <- NULL; CP <- NULL; EC <- NULL; SC <- NULL; AZ <- NULL
    for (i in 1:length(sample_names)) {
        GE <- datasetm[,i]
        Z1 <- (datasetm[match(IPSGdata$GENE, row.names(datasetm)),i]-mean(GE))/sd(GE)
        WEIGHT <- NULL; MIG <- NULL; k <- 1
        for (gen in unique(IPSGdata$NAME)) {
            MIG[k] <- mean(Z1[which(IPSGdata$NAME==gen)], na.rm=TRUE)
            WEIGHT[k] <- mean(IPSGdata$WEIGHT[which(IPSGdata$NAME==gen)])
            k<-k+1}
        WG <- MIG*WEIGHT
        MHC[i]<-mean(WG[1:10], na.rm = T)
        CP[i]<-mean(WG[11:20], na.rm = T)
        EC[i]<-mean(WG[21:24], na.rm = T)
        SC[i]<-mean(WG[25:26], na.rm = T)
        AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
        IPS[i]<-ipsmap(AZ[i])

        if(makeplot=="yes"){

            data_a <- data.frame (start = c(0, 2.5, 5, 7.5, 10, 15, seq(20, 39), 0, 10, 20, 30),
                                  end = c(2.5, 5, 7.5, 10, 15, seq(20, 40), 10, 20, 30, 40),
                                  y1 = c(rep(2.6, 26), rep(0.4, 4)), y2 = c(rep(5.6, 26), rep(2.2, 4)),
                                  z = c(MIG[c(21:26, 11:20, 1:10)], EC[i], SC[i], CP[i], MHC[i]),
                                  vcol = c(unlist(lapply(MIG[c(21:26, 11:20, 1:10)], mapcolors)),
                                           unlist(lapply(c(EC[i], SC[i], CP[i], MHC[i]),mapbw))),
                                  label = c(unique(IPSGdata$NAME)[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
            data_a$label <- factor(data_a$label, levels = unique(data_a$label))
            plot_a1 <- ggplot() +
                geom_rect(data_a, mapping = aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label),
                          size = 0.5, color = "black", alpha = 1) +
                coord_polar() + scale_y_continuous(limits = c(0, 6)) +
                scale_fill_manual(values = as.vector(data_a$vcol), guide = F) + theme_bw() +
                theme(panel.margin = unit(0, 'mm'), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), panel.border = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "white"),
                      axis.text = element_blank(), axis.ticks = element_blank()) +
                geom_text(aes(x = 5, y = 1.3, label = "EC"), size = 4) +
                geom_text(aes(x = 15, y = 1.3, label = "SC"), size = 4) +
                geom_text(aes(x = 25, y = 1.3, label = "CP"), size = 4) +
                geom_text(aes(x = 35, y = 1.3, label = "MHC"), size = 4)
            plot_a <- plot_a1 +
                geom_text(aes(x = c(1.25, 3.75, 6.25, 8.75, 17.5, 12.5), y = 4.1,
                              label = c("+ Act CD4", "+ Act CD8", "+ Tem CD4", "+ Tem CD8", "- MDSC", "- Treg")),
                          angle = c(78.75, 56.25, 33.75, 11.25, -67.5, -22.5), size = 4) +
                geom_text(aes(x = seq(20.5, 39.5, 1), y = 4.1,
                              label = c("PD-1 -","CTLA4 -","LAG3 -","TIGIT -","TIM3 -",
                                        "PD-L1 -","PD-L2 -","CD27 +","ICOS +","IDO1 -",
                                        "B2M +","TAP1 +","TAP2 +","HLA-A +","HLA-B +",
                                        "HLA-C +","HLA-DPA1 +","HLA-DPB1 +","HLA-E +","HLA-F +")),
                          angle = c(seq(85.5, 4.5, -9), seq(-4.5, -85.5, -9)), size = 4) +
                geom_text(aes(x = 0, y = 6, label = paste0("Immunophenoscore: ", IPS[i])),
                          angle = 0, size = 6, vjust = -0.5)+
                theme(axis.title = element_blank()) +
                theme(plot.margin = unit(c(0,0,0,0), "mm")) +
                geom_text(vjust = 1.15, hjust = 0, size = 4,
                          aes(x = 25.5, y = 6, hjust = 0,
                              label = "\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n"))

            legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"),
                                 panel.margin = unit(0,"null"), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),panel.border = element_blank(),
                                 panel.background = element_blank(), axis.text = element_blank(),
                                 axis.ticks = element_blank(), axis.title.x = element_blank(),
                                 axis.line = element_line(colour = "white"))
            data_bc <- data.frame(start = rep(0,23), end = rep(0.7,23), y1 = seq(0,22), y2 = seq(1,23),
                                  vcolb = c(unlist(lapply(seq(-3, 3, by = 6/22), mapcolors))),
                                  vcolc = c(unlist(lapply(seq(-2, 2, by = 4/22), mapbw))),
                                  label = LETTERS[1:23])
            ## Legend sample-wise (averaged) z-scores
            data_b_ticks <- data.frame(x = rep(1.2, 7), value = seq(-3, 3), y = seq(0, 6)*(22/6)+0.5)
            plot_b <- ggplot(hjust = 0) +
                geom_rect(data_bc, mapping = aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label),
                          size = 0.5, color = "black", alpha = 1) +
                scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) +
                scale_fill_manual(values = as.vector(data_bc$vcolb), guide = FALSE) +
                geom_text(data_b_ticks, mapping = aes(x = x, y = y, label = value), hjust = "inward", size = 4) +
                theme_bw() + legendtheme + ylab("Sample-wise (averaged) z-score")
            ## Legend weighted z-scores
            data_c_ticks <- data.frame(x = rep(1.2, 5), value = seq(-2, 2), y = seq(0, 4)*(22/4)+0.5)
            plot_c <- ggplot() +
                geom_rect(data_bc, mapping = aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label),
                          size = 0.5, color = "black", alpha = 1) +
                scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) +
                scale_fill_manual(values = as.vector(data_bc$vcolc), guide = FALSE) +
                geom_text(data_c_ticks, mapping = aes(x=x, y=y, label=value), hjust = "inward", size = 4) +
                theme_bw() + legendtheme + ylab("Weighted z-score")

            ## Save plot to file (1 pdf file for each sample)
            if(is.null(plotpath)){file_name <- paste0("IPS_", sample_names[i], ".pdf")
            } else {file_name <- paste0(plotpath, "/IPS_", sample_names[i], ".pdf")}
            pdf(file_name, width = 10, height = 8)
            grid.arrange(plot_a, plot_b, plot_c, ncol = 3, widths = c(0.8, 0.1, 0.1))
            dev.off()
        }
    }

    names(IPS) <- sample_names
    return(returnAsInput(userdata = dataset, result = IPS, SignName = "IPS"))
}


#' Core Matrisome Gene signature
#'
#' Given a dataset, it returns the median genes expression based on Yuzhalin et all. (2018).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples. Row names must be Official Symbol.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom matrixStats colMedians
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
MatriSign <- function(dataset, nametype) {

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        Matridata <- mapIds(org.Hs.eg.db, keys = Matridata, column = nametype, keytype = "SYMBOL", multiVals = "first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(Matridata %in% row.names(dataset)), " matrisome's genes out of 9\n"))

    median_cm <- colMedians(datasetm[row.names(datasetm) %in% Matridata, ])

    return(returnAsInput(userdata = dataset, result = median_cm, SignName = "MatrisomeScore"))
}

#' Mitotic Index
#'
#' Given a dataset, it returns the means genes expression based on Yang et all. (2016).
#'
#' @param dataset expression values where rows correspond to genes and columns correspond to samples. Row names must be Official Symbol.
#' @param nametype gene name ID of your dataset row names.
#'
#' @return NULL
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
MitoticIndexSign <- function(dataset, nametype) {

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!="SYMBOL"){
        MIdata <- mapIds(org.Hs.eg.db,keys= MIdata, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(MIdata %in% row.names(datasetm)), " mititotic index genes out of 9\n"))

    MI_means <- colMeans(datasetm[row.names(datasetm) %in% MIdata, ])

    return(returnAsInput(userdata = dataset, result = MI_means, SignName = "MitoticIndex"))
}
