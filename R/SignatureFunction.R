
matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)}

#' Endothelial-Mesenchymal Transition Signature
#'
#' Given a dataset, it returns the Endothelial score and the Mesenchymal score for each sample, based on QH Miow at all. (2015).
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

    Signature_EL <- EMTdata[grep('Epithelial-like', EMTdata$Category),]
    Signature_ML <- EMTdata[-grep('Epithelial-like', EMTdata$Category),]

    cat(paste0("The function is using ", sum(Signature_EL$Gene_Symbol %in% row.names(dataset)),
               " epithelial-like genes out of ", length(Signature_EL), "\nThe function is using ",
               sum(Signature_ML$Gene_Symbol %in% row.names(dataset)), " mesenchymal-like genes out of", length(Signature_ML),"\n"))

    gene_sets <- list(Epithelial=Signature_EL$Gene_Symbol, Mesenchimal=Signature_ML$Gene_Symbol)

    dots <- list(...)
    args <- matchArguments(dots, list(expr = as.matrix(dataset), gset.idx.list = gene_sets, method = "ssgsea",
                                      kcdf = "Poisson", abs.ranking = F, min.sz = 5, max.sz = Inf,
                                      parallel.sz = 1L, mx.diff = TRUE, ssgsea.norm = TRUE))
    gsva_matrix <- do.call(gsva, args)

    return(gsva_matrix)
}


#' Piroptosis Signature
#'
#' Given a dataset, it returns the piroptosis score for each sample, based on Mingjun Zheng et all. (2020).
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

    cat(paste0("The function is using ", sum(Pirodata$Gene_Symbol %in% row.names(dataset))," genes out of", length(Pirodata$Gene_Symbol)))
    Pirodata <- Pirodata[Pirodata$Gene_Symbol %in% row.names(dataset), ]
    Piroscore <- colSums(dataset[Pirodata$Gene_Symbol, ]*Pirodata$Coefficient)
    return(Piroscore)
}


#' FerroptosisSignature
#'
#' Given a dataset, it returns the Ferroptosis score for each sample Ying Ye et all. (2021).
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

    cat(paste0("The function is using ", sum(Ferrdata$Gene_Symbol %in% row.names(dataset))," genes out of", length(Ferrdata$Gene_Symbol)))
    Ferrdata <- Ferrdata[Ferrdata$Gene_Symbol %in% row.names(dataset), ]
    ferrscore <- colSums(dataset[Ferrdata$Gene_Symbol, ]*Ferrdata$Coefficient)
    return(ferrscore)
}


#' LIpidSignature
#'
#' Given a dataset, it returns the Lipid score for each sample Mingjun Zheng et all. (2020).
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

    cat(paste0("The function is using ", sum(Lipidata$Gene_Symbol %in% row.names(dataset))," genes out of", length(Lipidata$Gene_Symbol)))
    Lipidata <- Lipidata[Lipidata$Gene_Symbol %in% row.names(dataset), ]
    lipidscore <- colSums(dataset[Lipidata$Gene_Symbol, ] * Lipidata$Coefficient)
    return(lipidscore)
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

    cat(paste0("The function is using ", sum(genetouse %in% rownames(dataset)), " genes out of", length(Hypodata$Gene_Symbol)))
    dataset <- dataset[rownames(dataset) %in% genetouse, ]

    med_counts <- sort(setNames(colMedians(as.matrix(dataset)), colnames(dataset)))
    scores <- data.frame(E = med_counts, HS = scale(med_counts))

    return(scores)
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

    cat(paste("The function is using", sum(Platdata$up %in% row.names(dataset)),
              "up-genes out of", length(Platdata$up), "\nThe function is using", sum(Platdata$down %in% row.names(dataset)), "down-genes out of",length(Platdata$down),"\n"))
    dots <- list(...)
    args <- matchArguments(dots, list(expr = dataset, gset.idx.list = Platdata,
                                      method = "ssgsea", kcdf = "Poisson", min.sz=5))
    gsva_count <- do.call(gsva, args)

    return(gsva_count)
}


#' Prognostic Signature
#'
#' Given a dataset, it returns the Quantile assignation for each sample from J. Millstein et. all (2020).
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

    cat(paste("The function is using", sum(Progdata$Genes %in% row.names(dataset)),
              "genes out of", length(Progdata$Genes), "\n"))

    intergene <- intersect(row.names(dataset), names(Progdata$Genes))
    dataset <- dataset[intergene,]
    gene_coeff <- colSums(dataset*Progdata$Genes[intergene])

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

    return(quantile_prog)
}


#' Metabolic Signature
#'
#' Given a list of DEG, it returns a matrix with pathways score and a correspondent pvalue calculated with Bootstrapping. This signature is based on Rosario et. all (2018).
#'
#' @param DEdata matrix of differential expression genes where rows correspond to genes, first column correspond to Log2FoldChange and second column to its adjusted pvalue.
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
    return(cbind(path_score, pvals))
}
