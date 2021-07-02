
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

    signature <- read.csv("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/EMT/EMTsign.txt", header= T, sep = " ", row.names = 1)

    if(nametype!="SYMBOL"){
        signature$Gene_Symbol <- mapIds(org.Hs.eg.db,keys= signature$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    Signature_EL <- signature[grep('Epithelial-like', signature$Category),]
    Signature_ML <- signature[-grep('Epithelial-like', signature$Category),]

    cat(paste0("The function is using ", sum(Signature_EL$Gene_Symbol %in% row.names(dataset)),
               " epithelial-like genes\nThe function is using ",
               sum(Signature_ML$Gene_Symbol %in% row.names(dataset)), " mesenchymal-like genes\n"))

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

    pirosign <- read.csv("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/piroptosi/pirosign.txt", header= T, sep = " ")
    if(nametype!="SYMBOL"){
        pirosign$Gene_Symbol <- mapIds(org.Hs.eg.db,keys= pirosign$Gene_Symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    cat(paste0("The function is using ", sum(pirosign$Gene_Symbol %in% row.names(dataset))," genes"))
    pirosign <- pirosign[pirosign$Gene_Symbol %in% row.names(dataset), ]
    Piroscore <- colSums(dataset[pirosign$Gene_Symbol, ]*pirosign$Coefficient)
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

    fersign <- read.csv("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/Ferroptosis/Ferrsign.txt", header = T, sep = " ")
    if(nametype!="SYMBOL"){
        fersign$Gene_symbol <- mapIds(org.Hs.eg.db,keys= fersign$Gene_symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    cat(paste0("The function is using ", sum(fersign$Gene_symbol %in% row.names(dataset))," genes"))
    fersign <- fersign[fersign$Gene_symbol %in% row.names(dataset), ]
    ferrscore <- colSums(dataset[fersign$Gene_symbol, ]*fersign$score)
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

    lipsign <- read.csv("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/Lipid metabolism/LipSign.txt", sep = " ")
    if(nametype!="SYMBOL"){
        lipsign$Gene_symbol <- mapIds(org.Hs.eg.db,keys= lipsign$Gene_symbol, column= nametype, keytype="SYMBOL", multiVals="first")
    }

    cat(paste0("The function is using ", sum(lipsign$Gene_Symb %in% row.names(dataset))," genes"))
    lipsign <- lipsign[lipsign$Gene_Symb %in% row.names(dataset), ]
    lipidscore <- colSums(dataset[lipsign$Gene_Symb, ] * lipsign$Coeff)
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
#' @import org.Hs.eg.db
#'
#' @export
HypoSign <- function(dataset, nametype){

    if (!(nametype %in% c("SYMBOL","ENTREZID","ENSEMBL","ENSEMBLTRANS"))){
        stop("The name of genes must be either SYMBOL, ENTREZID, ENSEMBL or ENSEMBLTRANS")
    }

    Hyp_signature <- read.csv("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/ipossia/Hyp_signature.txt", sep = " ")

    if(nametype=="SYMBOL") { genetouse <- signature$genes
    } else if(nametype=="ENSEMBL") { genetouse <- signature$ensemble
    } else (Hyp_signature$Symbol <- mapIds(org.Hs.eg.db,keys= Hyp_signature$Symbol, column= nametype, keytype="SYMBOL", multiVals="first"))

    cat(paste0("The function is using ", sum(genetouse %in% rownames(dataset)), " genes out of 50"))
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

    load("/Users/Fabiola/Downloads/platinumres.RData")

    if(nametype!= "SYMBOL"){
        signature <- lapply(signature, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }


    cat(paste("The function is using", sum(signature$up %in% row.names(dataset)),
              "up-genes\nThe function is using", sum(signature$down %in% row.names(dataset)), "down-genes\n"))
    dots <- list(...)
    args <- matchArguments(dots, list(expr = dataset, gset.idx.list = signature,
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

    load("/Users/Fabiola/Downloads/platinumres.RData")

    if(nametype!="SYMBOL"){
        signature <- lapply(signature, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }

    intergene <- intersect(row.names(dataset), names(prognosticsign$Genes))
    dataset <- dataset[intergene,]
    gene_coeff <- colSums(dataset*prognosticsign$Genes[intergene])

    age_coef <- sapply(age, function(p)
        if(p<=53){0} else if(p>53 & p<=60){prognosticsign$Age[1]
        } else if(p>60 & p<=67){prognosticsign$Age[2]} else {prognosticsign$Age[3]})

    stage_coef <- sapply(stage, function(p)
        if(p=="NA"){prognosticsign$Stage[2]} else if(p=="I"|p=="II"){prognosticsign$Stage[1]} else {0})

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

    load("/Users/Fabiola/Desktop/TESI/R file/VIOlets analysis/metabolic dysregulation/MetabolicPathways.RData")

    if(nametype!="SYMBOL"){
        metabolicpath <- lapply(metabolicpath, function(x)
            suppressMessages(mapIds(org.Hs.eg.db,keys= x, column= nametype, keytype="SYMBOL", multiVals="first")))
    }

    gene_pathway <- lapply(metabolicpath, intersect, row.names(DEdata))
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
