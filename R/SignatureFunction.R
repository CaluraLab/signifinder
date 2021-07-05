
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
            if(is.vector(result)){userdata@meta.data <- cbind(userdata@meta.data, SignName=result)
            } else {userdata@meta.data <- cbind(userdata@meta.data, t(result))}
        } else if(class(userdata)%in%c("SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment")){
            if(is.vector(result)){userdata@colData <- cbind(userdata@colData, SignName=result)
            } else {userdata@colData <- cbind(userdata@colData, t(result))}}
        return(userdata)
    } else {return(result)}}

#' Endothelial-Mesenchymal Transition Signature
#'
#' Given a dataset, it returns the Endothelial score and the Mesenchymal score for each sample, based on QH Miow at al. (2015).
#'
#' @param dataset a matrix of expression values where rows correspond to genes and columns correspond to samples. Row names must be Official Symbol.
#' @param nametype gene name ID of your dataset row names.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
EMTSign <- function(dataset, nametype, ...) {

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
               " epithelial-like genes out of ", length(Signature_EL), "\nThe function is using ",
               sum(Signature_ML$Gene_Symbol %in% row.names(datasetm)), " mesenchymal-like genes out of",
               length(Signature_ML),"\n"))

    gene_sets <- list(Epithelial=Signature_EL$Gene_Symbol, Mesenchimal=Signature_ML$Gene_Symbol)

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = gene_sets, method = "ssgsea",
                                      kcdf = "Poisson", abs.ranking = F, min.sz = 5, max.sz = Inf,
                                      parallel.sz = 1L, mx.diff = TRUE, ssgsea.norm = TRUE))
    gsva_matrix <- do.call(gsva, args)

    return(returnAsInput(userdata = dataset, result = gsva_matrix, SignName = ""))
}


#' Piroptosis Signature
#'
#' Given a dataset, it returns the piroptosis score for each sample, based on Mingjun Zheng et al. (2020).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples. Row names must be Official Symbol.
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
        Pirodata$Gene_Symbol <- mapIds(org.Hs.eg.db,keys= Pirodata$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(Pirodata$Gene_Symbol %in% row.names(datasetm))," genes out of", length(Pirodata$Gene_Symbol)))
    Pirodata <- Pirodata[Pirodata$Gene_Symbol %in% row.names(datasetm), ]
    Piroscore <- colSums(datasetm[Pirodata$Gene_Symbol, ]*Pirodata$Coefficient)
    return(returnAsInput(userdata = dataset, result = Piroscore, SignName = "PiroptosisScore"))
}


#' FerroptosisSignature
#'
#' Given a dataset, it returns the Ferroptosis score for each sample Ying Ye et al. (2021).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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

    cat(paste0("The function is using ", sum(Ferrdata$Gene_Symbol %in% row.names(datasetm))," genes out of", length(Ferrdata$Gene_Symbol)))
    Ferrdata <- Ferrdata[Ferrdata$Gene_Symbol %in% row.names(datasetm), ]
    ferrscore <- colSums(datasetm[Ferrdata$Gene_Symbol, ]*Ferrdata$Coefficient)
    return(returnAsInput(userdata = dataset, result = ferrscore, SignName = "FerroptosisScore"))
}


#' LIpidSignature
#'
#' Given a dataset, it returns the Lipid score for each sample Mingjun Zheng et al. (2020).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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

    cat(paste0("The function is using ", sum(Lipidata$Gene_Symbol %in% row.names(datasetm))," genes out of", length(Lipidata$Gene_Symbol)))
    Lipidata <- Lipidata[Lipidata$Gene_Symbol %in% row.names(datasetm), ]
    lipidscore <- colSums(datasetm[Lipidata$Gene_Symbol, ] * Lipidata$Coefficient)
    return(returnAsInput(userdata = dataset, result = lipidscore, SignName = "LipidScore"))
}


#' Hypoxia Signature
#'
#' Given a dataset, it returns the hypoxia score for each sample as in Buffa et al. 2010.
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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
    } else (genetouse <- mapIds(org.Hs.eg.db,keys= Hypodata$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first"))

    datasetm <- getMatrix(dataset)

    cat(paste0("The function is using ", sum(genetouse %in% rownames(datasetm)),
               " genes out of", length(Hypodata$Gene_Symbol)))
    datasetm <- datasetm[rownames(datasetm) %in% genetouse, ]

    med_counts <- sort(setNames(colMedians(as.matrix(datasetm)), colnames(datasetm)))
    scores <- data.frame(E = med_counts, HS = scale(med_counts))

    return(returnAsInput(userdata = dataset, result = scores, SignName = ""))
}


#' Platinum Resistance Signature
#'
#' Given a dataset, it returns the gsva score for each sample from International Cancer Genome Consortium (ICGC).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param nametype gene name ID of your dataset row names.
#' @param ... other arguments passed on to the GSVA function.
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @export
PlatResSign <- function(dataset, nametype,  ...){

    if (!nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS")){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    if(nametype!= "SYMBOL"){
        Platdata <- lapply(Platdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(Platdata$up %in% row.names(datasetm)),
              "up-genes out of", length(Platdata$up), "\nThe function is using",
              sum(Platdata$down %in% row.names(datasetm)), "down-genes out of", length(Platdata$down),"\n"))

    dots <- list(...)
    args <- matchArguments(dots, list(expr = datasetm, gset.idx.list = Platdata,
                                      method = "ssgsea", kcdf = "Poisson", min.sz=5))
    gsva_count <- do.call(gsva, args)

    return(returnAsInput(userdata = dataset, result = gsva_count, SignName = ""))
}


#' Prognostic Signature
#'
#' Given a dataset, it returns the Quantile assignation for each sample from J. Millstein et. al (2020).
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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
        Progdata <- lapply(Progdata, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }

    datasetm <- getMatrix(dataset)

    cat(paste("The function is using", sum(Progdata$Genes %in% row.names(datasetm)),
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

    return(returnAsInput(userdata = dataset, result = quantile_prog, SignName = "PrognosticScore"))
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
    return(returnAsInput(userdata = dataset, result = rbind(path_score, pvals), SignName = ""))
}


#' Immunogenic Signature
#'
#' Given a dataset, it returns the ImmunoScore for each sample. This signature is based on Dapeng Hao et. al (2018).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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
        ImmunoGenes$genes <- mapIds(org.Hs.eg.db, keys = ImmunoGenes$genes, column = nametype,
                                    keytype = "SYMBOL", multiVals = "first")
    }

    datasetm <- getMatrix(dataset)

    g <- intersect(row.names(datasetm), ImmunoGenes$genes)

    subdataset <- datasetm[g,]
    ImmunoGenes <- ImmunoGenes[ImmunoGenes$genes %in% g, ]

    SE <- (ImmunoGenes$HR - ImmunoGenes$`95CI_L`)/1.96
    k <- (1 - ImmunoGenes$HR)/SE

    ImmunoScores <- unlist(lapply(seq_len(ncol(subdataset)), function(p) sum(k*subdataset[,p], na.rm = T)))

    return(returnAsInput(userdata = dataset, result = ImmunoScores, SignName = "ImmunoScore"))
}


#' ConsensusOV Signature
#'
#' Given a dataset, it returns ovarian cancer subtypes. This signature is based on Chen et. al (2018).
#'
#' @param dataset matrix of expression values where rows correspond to genes and columns correspond to samples.
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
