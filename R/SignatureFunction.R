

#' Epithelial-Mesenchymal Transition Signature
#'
#' @description This signature is computed accordingly to the reference paper,
#' to have more details explore the function
#' \code{\link[signifinder]{availableSignatures}}.
#'
#' @param dataset Normalized expression values. A data frame or a matrix
#' where rows correspond to genes and columns correspond to samples.
#' Alternatively, an object of type \linkS4class{SummarizedExperiment},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}} where the normalized
#' expression values should be in an assay called 'norm_expr'.
#' @param nametype character string saying the type of gene name ID (row names
#' in dataset). Either one of "SYMBOL", "ENTREZID" or "ENSEMBL".
#' @param inputType character string saying the type of data you are using.
#' Either one of "microarray" or "rnaseq".
#' @param author character string saying the first author of the signature
#' publication. Check it in \code{\link[signifinder]{availableSignatures}}.
#' @param whichAssay integer scalar or string indicating which assay of
#' dataset to use.
#' @param hgReference character string saying the human reference genome.
#' Either one of "hg19" or "hg38".
#' @param isMalignant logical vector of the same length of ncol(dataset), where
#' TRUE states malignant cells and FALSE states non-malignant cells.
#' @param ... other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
#' function.
#'
#' @return If dataset is a \linkS4class{SummarizedExperiment} object, then
#' scores are added in the \code{\link[SummarizedExperiment]{colData}} section.
#' If dataset is a data frame or a matrix, then a
#' \linkS4class{SummarizedExperiment} object is created in which scores are
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom stats prcomp
#'
#' @examples
#' data(ovse)
#' EMTSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
EMTSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    author = "Miow", whichAssay = "norm_expr", hgReference = "hg38",
    isMalignant = NULL, ...) {

    .consistencyCheck(nametype, "EMTSign", author)
    datasetm <- .getMatrix(dataset, whichAssay)

    if (author == "Miow") {
        sign_df <- EMT_Miow
        sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

        EL <- sign_df[grep("Epithelial-like", sign_df$class), ]
        ML <- sign_df[grep("Mesenchymal-like", sign_df$class), ]

        .percentageOfGenesUsed(
            "EMTSign", datasetm, EL$SYMBOL, "epithelial", author = author)
        .percentageOfGenesUsed(
            "EMTSign", datasetm, ML$SYMBOL, "mesenchymal", author = author)

        gene_sets <- list( Epithelial = EL$SYMBOL, Mesenchymal = ML$SYMBOL)
        names(gene_sets) <- paste0("EMT_Miow_", names(gene_sets))

        dots <- list(...)
        args <- .matchArguments(dots, list(
            exprData = datasetm, geneSets = gene_sets, normalize = FALSE))
        gsvaPar <- do.call(ssgseaParam, args)
        score <- gsva(gsvaPar, verbose = FALSE)

        return(.returnAsInput(
            userdata = dataset, result = score, SignName = "", datasetm))
    } else if (author == "Barkley") {

        .isMalignantCheck(isMalignant, dataset)

        sign_df <- PanState_Barkley
        sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

        .percentageOfGenesUsed(
            "EMTSign", datasetm,
            sign_df$SYMBOL[sign_df$class == "cEMT"], "cEMT")
        .percentageOfGenesUsed(
            "EMTSign", datasetm,
            sign_df$SYMBOL[sign_df$class == "pEMT"], "pEMT")

        sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
        sign_list <- split(sign_df$SYMBOL, sign_df$class)
        names(sign_list) <- paste0("EMT_Barkley_", names(sign_list))

        datasetm_n <- datasetm[,isMalignant]
        modules = c("EMT_Barkley_cEMT", "EMT_Barkley_pEMT")
        mod <- sign_list[modules]
        score <- .barkleyFun(
            dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)
        colnames(score) <- colnames(datasetm)

        return(.returnAsInput(
            userdata = dataset, result = score, SignName = "", datasetm))
    } else {
        if (author == "Mak") {
            sign_df <- EMT_Mak
            sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

            Sign_E <- sign_df$SYMBOL[sign_df$class == "E"]
            Sign_M <- sign_df$SYMBOL[sign_df$class == "M"]

            .percentageOfGenesUsed(
                "EMTSign", datasetm, Sign_E, "epithelial", author = author)
            .percentageOfGenesUsed(
                "EMTSign", datasetm, Sign_M, "mesenchymal", author = author)

            columnNA <- .managena(datasetm, genes = sign_df$SYMBOL)
            score <- colMeans(
                datasetm[intersect(Sign_M, row.names(datasetm)), ]) - colMeans(
                    datasetm[intersect(Sign_E, row.names(datasetm)), ])
            score[columnNA > 0.9] <- NA
        } else if (author == "Cheng") {
            sign_df <- EMT_Cheng
            sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

            .percentageOfGenesUsed(
                "EMTSign", datasetm, sign_df$SYMBOL, author = author)

            datasetm_n <- datasetm[
                intersect(row.names(datasetm), sign_df$SYMBOL), ]

            datasetm_n <- if (inputType == "rnaseq") {
                log2(datasetm_n + 1)
            } else { datasetm_n }

            columnNA <- .managena(datasetm = datasetm_n, genes = sign_df$SYMBOL)

            score <- prcomp(t(datasetm_n))$x[, 1]
            score[columnNA > 0.9] <- NA
        } else if (author == "Thompson") {
            sign_df <- EMT_Thompson
            sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

            Sign_E <- sign_df$SYMBOL[sign_df$class == "epithelial"]
            Sign_M <- sign_df$SYMBOL[sign_df$class == "mesenchymal"]

            .percentageOfGenesUsed(
                "EMTSign", datasetm, Sign_E, "epithelial", author = author)
            .percentageOfGenesUsed(
                "EMTSign", datasetm, Sign_M, "mesenchymal", author = author)

            dataset <- .dataTransformation(
                dataset = dataset, data = datasetm, trans = "CPM",
                hgReference = hgReference, nametype = nametype)
            datasetm_n <- as.matrix(assays(dataset)[["CPM"]])
            datasetm_n <- log2(datasetm_n + 1)

            t_dataset <- t(datasetm_n)
            epi <- scale(t_dataset[ ,intersect(Sign_E, colnames(t_dataset))])
            mes <- scale(t_dataset[ ,intersect(Sign_M, colnames(t_dataset))])
            epi <- rowSums(log2(epi - min(epi) + 1))
            mes <- rowSums(log2(mes - min(mes) + 1))
            score <- mes-epi
        }
        return(.returnAsInput(
            userdata = dataset, result = score,
            SignName = paste0("EMT_", author), datasetm))
    }
}


#' Pyroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' pyroptosisSign(dataset = ovse)
#'
#' @export
pyroptosisSign <- function(
    dataset, nametype = "SYMBOL", inputType = "rnaseq",
    author = "Ye", whichAssay = "norm_expr", hgReference = "hg38") {

    .consistencyCheck(nametype, "pyroptosisSign", author)

    sign_df <- get(paste0("Pyroptosis_", author))
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- .getMatrix(dataset, whichAssay)

    if (author == "Ye") {
        dataset <- .dataTransformation(
            dataset, datasetm, "FPKM", hgReference, nametype)
        datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
        datasetm_n <- scale(datasetm_n)
    } else if (author == "Shao") {
        if (inputType == "rnaseq") {
            dataset <- .dataTransformation(
                dataset, datasetm, "FPKM", hgReference, nametype)
            datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
        } else { datasetm_n <- datasetm }
    } else if (author == "Lin") {
        dataset <- .dataTransformation(
            dataset, datasetm, "TPM", hgReference, nametype)
        datasetm_n <- as.matrix(assays(dataset)[["TPM"]])
    } else { datasetm_n <- datasetm}

    score <- .coeffScore(sign_df, datasetm_n, "pyroptosisSign", author = author)

    return(.returnAsInput(
        userdata = dataset, result = score,
        SignName = paste0("Pyroptosis_", author), datasetm))
}


#' Ferroptosis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' ferroptosisSign(dataset = ovse)
#'
#' @export
ferroptosisSign <- function(
    dataset, nametype = "SYMBOL", inputType = "rnaseq",
    author = "Ye", whichAssay = "norm_expr", hgReference = "hg38") {

  .consistencyCheck(nametype, "ferroptosisSign", author)

  sign_df <- get(paste0("Ferroptosis_", author))
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  if (author %in% c("Liu", "Ye") & inputType == "rnaseq") {
    dataset <- .dataTransformation(
      dataset, datasetm, "FPKM", hgReference, nametype)
    datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
    datasetm_n <- datasetm_n[rownames(datasetm_n) %in% sign_df$SYMBOL, ]
  } else {
    datasetm_n <- datasetm[rownames(datasetm) %in% sign_df$SYMBOL, ]
  }

  score <- .coeffScore(
    sign_df, datasetm_n, "ferroptosisSign", author = author)

  if (author == "Liang") { score <- exp(score) }

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("Ferroptosis_", author), datasetm))
}


#' Lipid Metabolism Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' lipidMetabolismSign(dataset = ovse)
#'
#' @export
lipidMetabolismSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "lipidMetabolismSign")

  sign_df <- LipidMetabolism_Zheng
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)
  score <- .coeffScore(sign_df, datasetm, "lipidMetabolismSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "LipidMetabolism_Zheng", datasetm))
}


#' Hypoxia Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' hypoxiaSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
hypoxiaSign <- function(
    dataset, nametype = "SYMBOL", author =  "Buffa", inputType = "microarray",
    whichAssay = "norm_expr", isMalignant = NULL) {

  .consistencyCheck(nametype, "hypoxiaSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Buffa") {
    datasetm_n <- if (inputType == "rnaseq") {
      log2(datasetm + 1)
    } else { datasetm }
    score <- .statScore(
      Hypoxia_Buffa$SYMBOL, datasetm = abs(datasetm_n), nametype = nametype,
      typeofstat = "median", namesignature = "hypoxiaSign")
    return(.returnAsInput(
      userdata = dataset, result = as.vector(scale(score)),
      SignName = "Hypoxia_Buffa", datasetm))
  } else {
    if (author == "Barkley") {

      .isMalignantCheck(isMalignant, dataset)

      sign_df <- PanState_Barkley
      sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

      .percentageOfGenesUsed(
        "hipoxiaSign", datasetm,
        sign_df$SYMBOL[sign_df$class == "Hypoxia"], "Hypoxia")

      sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
      sign_list <- split(sign_df$SYMBOL, sign_df$class)
      names(sign_list) <- paste0(names(sign_list), "_Barkley")

      datasetm_n <- datasetm[,isMalignant]
      modules = "Hypoxia_Barkley"
      mod <- sign_list[modules]
      score <- .barkleyFun(
          dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

      return(.returnAsInput(
        userdata = dataset, result = score,
        SignName = "Hypoxia_Barkley", datasetm))}
  }
}


#' Immunogenic Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' immunoScoreSign(dataset = ovse)
#'
#' @export
immunoScoreSign <- function(
    dataset, nametype = "SYMBOL", author = "Hao", inputType = "rnaseq",
    whichAssay = "norm_expr", hgReference = "hg38") {

  .consistencyCheck(nametype, "immunoScoreSign", author)
  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Hao") {
    sign_df <- ImmunoScore_Hao
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    g <- intersect(row.names(datasetm), sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "immunoScoreSign", datasetm, sign_df$SYMBOL, author = author)

    if (inputType == "rnaseq") {
      dataset <- .dataTransformation(
        dataset, datasetm, "FPKM", hgReference, nametype)
      datasetm_n <- log2(as.matrix(assays(dataset)[["FPKM"]]) + 0.01)
      datasetm_n <- datasetm_n[g, ]
    } else { datasetm_n <- datasetm[g, ] }

    sign_df <- sign_df[sign_df$SYMBOL %in% g, ]
    columnNA <- .managena(datasetm_n, g)
    SE <- (sign_df$HR - sign_df$`95CI_L`) / 1.96
    k <- (1 - sign_df$HR) / SE
    score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(p) {
      sum(k * datasetm_n[, p], na.rm = TRUE) }))
    score[columnNA > 0.9] <- NA
  } else if (author == "Roh") {
    datasetm_n <- log2(datasetm + 1)
    score <- .statScore(
      ImmunoScore_Roh$SYMBOL, datasetm = datasetm_n,
      nametype = nametype, typeofstat = ".meang",
      namesignature = "immunoScoreSign", author = author)
  }

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("ImmunoScore_", author), datasetm))
}


#' ConsensusOV Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... optional parameters to be passed to
#' \code{\link[consensusOV]{get.subtypes}}.
#'
#' @inherit EMTSign return
#'
#' @importFrom consensusOV get.subtypes
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @examples
#' data(ovse)
#' consensusOVSign(dataset = ovse)
#'
#' @export
consensusOVSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr", ...) {

  .consistencyCheck(nametype, "consensusOVSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  if (nametype != "ENTREZID") {
    genename <- mapIds(
      org.Hs.eg.db,
      keys = row.names(datasetm),
      column = "ENTREZID",
      keytype = nametype,
      multiVals = "first")
    datasetm_n <- datasetm[!is.na(genename), ]
    genename <- genename[!is.na(genename)]
    datasetm_n <- datasetm_n[!duplicated(genename), ]
    genename <- genename[!duplicated(genename)]
  } else {
    genename <- row.names(datasetm)
    datasetm_n <- datasetm
  }

  consensus_subtypes <- get.subtypes(
    expression.dataset = datasetm_n, entrez.ids = genename, ...)
  scores <- consensus_subtypes$rf.probs
  colnames(scores) <- paste0(
    "ConsensusOV_Chen_", substring(colnames(scores), 1, 3))

  return(.returnAsInput(
    userdata = dataset, result = t(scores), SignName = "", datasetm))
}


#' ImmunoPhenoScore Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom stats sd
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' IPSSign(dataset = ovse)
#'
#' @export
IPSSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    hgReference = "hg38") {
  ## code adapted from https://github.com/icbi-lab/Immunophenogram

  .consistencyCheck(nametype, "IPSSign")

  sign_df <- IPS_Charoentong
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)
  sample_names <- colnames(datasetm)

  .percentageOfGenesUsed("IPSSign", datasetm, sign_df$SYMBOL)

  dataset <- .dataTransformation(
    dataset, datasetm, "TPM", hgReference, nametype)
  datasetm_n <- log2(as.matrix(assays(dataset)[["TPM"]]) + 1)
  datasetm_n <- datasetm_n[rownames(datasetm_n) %in% sign_df$SYMBOL, ]

  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL
  for (i in seq_len(length(sample_names))) {
    GE <- datasetm_n[, i]
    Z1 <- (datasetm_n[match(sign_df$SYMBOL, row.names(datasetm_n)), i] -
             mean(GE)) / sd(GE)
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    for (gen in unique(sign_df$NAME)) {
      MIG[k] <- mean(Z1[which(sign_df$NAME == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(sign_df$WEIGHT[which(sign_df$NAME == gen)])
      k <- k + 1}

    WG <- MIG * WEIGHT
    MHC[i] <- mean(WG[which(
      unique(sign_df$NAME) %in% unique(
        sign_df$NAME[sign_df$class == "MHC"]))], na.rm = TRUE)
    CP[i] <- mean(WG[which(
      unique(sign_df$NAME) %in% unique(
        sign_df$NAME[sign_df$class == "CP"]))], na.rm = TRUE)
    EC[i] <- mean(WG[which(
      unique(sign_df$NAME) %in% unique(
        sign_df$NAME[sign_df$class == "EC"]))], na.rm = TRUE)
    SC[i] <- mean(WG[which(
      unique(sign_df$NAME) %in% unique(
        sign_df$NAME[sign_df$class == "SC"]))], na.rm = TRUE)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- .ipsmap(AZ[i])
  }

  ipsres <- data.frame(IPS, MHC, CP, EC, SC)
  rownames(ipsres) <- sample_names
  colnames(ipsres) <- c(
    "IPS_Charoentong", "IPS_Charoentong_MHC", "IPS_Charoentong_CP",
    "IPS_Charoentong_EC", "IPS_Charoentong_SC")

  return(.returnAsInput(
    userdata = dataset, result = t(ipsres), SignName = "", datasetm))
}


#' Core Matrisome Gene signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' matrisomeSign(dataset = ovse)
#'
#' @export
matrisomeSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "matrisomeSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  score <- .statScore(
    Matrisome_Yuzhalin$SYMBOL, datasetm,
    nametype, "median", "matrisomeSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Matrisome_Yuzhalin", datasetm))
}


#' Mitotic Index
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' mitoticIndexSign(dataset = ovse)
#'
#' @export
mitoticIndexSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "mitoticIndexSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  score <- .statScore(
    MitoticIndex_Yang$SYMBOL, datasetm,
    nametype, "mean", "mitoticIndexSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "MitoticIndex_Yang", datasetm))
}


#' Immune Cytolytic Activity Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' immuneCytSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
immuneCytSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    author = "Rooney", whichAssay = "norm_expr", hgReference = "hg38") {

  .consistencyCheck(nametype, "immuneCytSign", author)
  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Rooney") {
    if (inputType == "rnaseq") {
      dataset <- .dataTransformation(
        dataset, datasetm, "TPM", hgReference, nametype)
      datasetm_n <- as.matrix(assays(dataset)[["TPM"]]) + 0.01
    } else { datasetm_n <- datasetm }
    score <- .statScore(
      ImmuneCyt_Rooney$SYMBOL, datasetm_n, nametype,
      ".meang", "ImmuneCytSign", author = author)
  } else if (author == "Davoli") {
    datasetm_n <- datasetm[
      row.names(datasetm) %in% ImmuneCyt_Davoli$SYMBOL, ]
    datasetm_r <- apply(datasetm_n, 1, rank)
    datasetm_r <- (datasetm_r - 1) / (nrow(datasetm_r) - 1)
    score <- .statScore(
      ImmuneCyt_Davoli$SYMBOL, t(datasetm_r), nametype,
      ".meang", "ImmuneCytSign", author = author)
  }
  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("ImmuneCyt_", author), datasetm))
}


#' IFN-gamma Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' IFNSign(dataset = ovse)
#'
#' @export
IFNSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "IFNSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  datasetm_n <- log10(datasetm + 1)
  score <- .statScore(
    IFN_Ayers$SYMBOL, datasetm_n, nametype, "mean", "IFNSign")

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "IFN_Ayers", datasetm))
}


#' ExpandedImmune Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' expandedImmuneSign(dataset = ovse)
#'
#' @export
expandedImmuneSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "expandedImmuneSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  datasetm_n <- log10(datasetm + 1)
  score <- .statScore(
    ExpandedImmune_Ayers$SYMBOL, datasetm_n,
    nametype, "mean", "expandedImmuneSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "ExpandedImmune_Ayers", datasetm))
}


#' TinflamSign Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' TinflamSign(dataset = ovse)
#'
#' @export
TinflamSign <- function(
    dataset, nametype = "SYMBOL", author= "Ayers", whichAssay = "norm_expr",
    hgReference = "hg38") {

  .consistencyCheck(nametype, "TinflamSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Ayers"){
    sign_df <- Tinflam_Ayers
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm_n <- log10(datasetm + 1)

    housekeeping <- intersect(
      row.names(datasetm_n), sign_df$SYMBOL[sign_df$class == "Housekeeping"])
    genes_pred <- intersect(
      row.names(datasetm_n), sign_df$SYMBOL[sign_df$class == "TInflam"])

    housekeeping_m <- apply(datasetm_n[housekeeping, ], 2, mean)
    datasetm_n <- sweep(datasetm_n[genes_pred, ], 2, housekeeping_m, FUN = "-")
    score <- .coeffScore(
      sign_df[sign_df$SYMBOL %in% genes_pred, ], datasetm_n, "TinflamSign")}

  else if (author == "Thompson"){
    sign_df <- Tinflam_Thompson
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "TinflamSign", datasetm, sign_df$SYMBOL, author = author)

    dataset <- .dataTransformation(
      dataset = dataset, data = datasetm, trans = "CPM",
      hgReference = hgReference, nametype = nametype)
    datasetm_n <- as.matrix(assays(dataset)[["CPM"]])
    datasetm_n <- log2(datasetm_n + 1)

    t_dataset <- t(datasetm_n)
    sc_dataset <- scale(
      t_dataset[, intersect(sign_df$SYMBOL, colnames(t_dataset))])
    score <- rowSums(log2(sc_dataset - min(sc_dataset) + 1))}

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("Tinflam_", author), datasetm))
}


#' Tertiary Lymphoid Structures (TLS) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' TLSSign(dataset = ovse)
#'
#' @export
TLSSign <- function(
    dataset, nametype = "SYMBOL", inputType = "rnaseq",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "TLSSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- if (inputType == "rnaseq") {
    log2(datasetm + 1) } else { datasetm }
  score <- .statScore(
    TLS_Cabrita$SYMBOL, datasetm_n, nametype, "mean", "TLSSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "TLS_Cabrita", datasetm))
}


#' CD49fHi Basal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' stemCellCD49fSign(dataset = ovse)
#'
#' @export
stemCellCD49fSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {
  .consistencyCheck(nametype, "StemCellCD49fSign")

  sign_df <- StemCellCD49f_Smith
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- log2(datasetm + 1)

  score <- .coeffScore(sign_df, datasetm_n, "StemCellCD49fSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "StemCellCD49f_Smith", datasetm))
}


#' Chromosomal instability Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' CINSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
CINSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "CINSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- if (inputType == "rnaseq") {
    log2(datasetm + 1)
  } else { datasetm }
  datasetm_n <- scale(datasetm_n, center = TRUE, scale = FALSE)

  score25 <- .statScore(
    CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"],
    datasetm_n, nametype, "sum", "CINSign")

  score70 <- .statScore(
    CIN_Carter$SYMBOL, datasetm_n, nametype, "sum", "CINSign")

  score <- data.frame(score25, score70)
  colnames(score) <- c("CIN_Carter_25", "CIN_Carter_70")

  return(.returnAsInput(
    userdata = dataset, result = t(score), SignName = "", datasetm))
}


#' Cell-cycle Signature classifier
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' cellCycleSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
cellCycleSign <- function(
    dataset, nametype = "SYMBOL", author = "Lundberg",
    inputType = "microarray", whichAssay = "norm_expr", isMalignant = NULL) {

  .consistencyCheck(nametype, "cellCycleSign", author)
  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Lundberg") {
    datasetm_n <- if (inputType == "rnaseq") {
      log2(datasetm + 1)
    } else { datasetm }
    score <- .statScore(
      CellCycle_Lundberg$SYMBOL, datasetm_n, nametype,
      "sum", "cellCycleSign", author = author)
  } else if (author == "Davoli") {
    datasetm_n <- datasetm[
      row.names(datasetm) %in% CellCycle_Davoli$SYMBOL, ]
    datasetm_r <- apply(datasetm_n, 1, rank)
    datasetm_r <- (datasetm_r - 1) / (nrow(datasetm_r) - 1)
    score <- .statScore(
      CellCycle_Davoli$SYMBOL, t(datasetm_r), nametype,
      ".meang", "cellCycleSign", author = author)
  } else if (author == "Barkley") {

    .isMalignantCheck(isMalignant, dataset)

    sign_df <- PanState_Barkley
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "cellCycleSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Cycle"], "Cycle")

    sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
    sign_list <- split(sign_df$SYMBOL, sign_df$class)
    names(sign_list) <- paste0("Cell", names(sign_list), "_Barkley")

    datasetm_n <- datasetm[,isMalignant]
    modules = "CellCycle_Barkley"
    mod <- sign_list[modules]
    score <- .barkleyFun(
        dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)
  }
  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("CellCycle_", author), datasetm))
}


#' Chemokine Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom stats prcomp
#'
#' @examples
#' data(ovse)
#' chemokineSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
chemokineSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "chemokineSign")

  sign_df <- Chemokines_Messina
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  .percentageOfGenesUsed("ChemokineSign", datasetm, sign_df$SYMBOL)

  datasetm_n <- datasetm[intersect(row.names(datasetm), sign_df$SYMBOL), ]
  datasetm_n <- if (inputType == "rnaseq") {
    log2(datasetm_n + 1)
  } else { datasetm_n }
  columnNA <- .managena(datasetm_n, sign_df$SYMBOL)
  score <- prcomp(t(datasetm_n), center = TRUE, scale = TRUE)$x[, 1]
  score[columnNA > 0.9] <- NA

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Chemokines_Messina", datasetm))
}


#' Adult Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' ASCSign(dataset = ovse)
#'
#' @export
ASCSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {
  .consistencyCheck(nametype, "ASCSign")

  sign_df <- ASC_Smith
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  .percentageOfGenesUsed("ASCSign", datasetm, sign_df$SYMBOL)

  datasetm_n <- log2(datasetm[row.names(datasetm) %in% sign_df$SYMBOL, ] + 1)

  columnNA <- .managena(datasetm_n, sign_df$SYMBOL)

  score <- rowSums(
    scale(t(datasetm_n), center = TRUE, scale = TRUE), na.rm = TRUE)
  score[columnNA > 0.9] <- NA

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "ASC_Smith", datasetm))
}


#' passON Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom SummarizedExperiment assays
#' @importFrom stats weighted.mean
#'
#' @examples
#' data(ovse)
#' PassONSign(dataset = ovse)
#'
#' @export
PassONSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    hgReference = "hg38", ...) {

  .consistencyCheck(nametype, "PassONSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  dataset <- .dataTransformation(
    dataset, datasetm, "TPM", hgReference, nametype)
  datasetm_n <- as.matrix(assays(dataset)[["TPM"]])

  sign_df <- PassON_Du
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)

  .percentageOfGenesUsed("PassONSign", datasetm_n, sign_df$SYMBOL)

  dots <- list(...)
  args <- .matchArguments(dots, list(
    exprData = datasetm_n, geneSets = sign_list))

  gsvaPar <- do.call(ssgseaParam, args)
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)

  gsva_mean <- vapply(seq_len(ncol(gsva_matrix)), function(x) {
    weighted.mean(gsva_matrix[, x], lengths(sign_list))},
    double(1))

  return(.returnAsInput(
    userdata = dataset, result = gsva_mean,
    SignName = "PassON_Du", datasetm))
}


#' IPRES Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA gsva ssgseaParam
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' IPRESSign(dataset = ovse)
#'
#' @export
IPRESSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    hgReference = "hg38", ...) {

  .consistencyCheck(nametype, "IPRESSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  dataset <- .dataTransformation(
    dataset, datasetm, "CPM", hgReference, nametype)
  datasetm_n <- log2(as.matrix(assays(dataset)[["CPM"]]))

  sign_df <- IPRES_Hugo
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)

  .percentageOfGenesUsed("IPRESSign", datasetm_n, sign_df$SYMBOL)

  dots <- list(...)
  args <- .matchArguments(dots, list(
    exprData = datasetm_n, geneSets = sign_list))

  gsvaPar <- do.call(ssgseaParam, args)
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
  score <- rowMeans(vapply(
    as.data.frame(t(gsva_matrix)), scale, double(ncol(gsva_matrix))))

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "IPRES_Hugo", datasetm))
}


#'  CIS (carcinoma-in situ) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' data(ovse)
#' CISSign(dataset = ovse)
#'
#' @export
CISSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {
  ## code instruction from an e-mail by Jaegil Kim:
  ## 'We first compute "log2(fold change) = log2(expression) - log2(median
  ## expression across samples)" for each gene and define the signature
  ## score as "mean of log2(fold change) across CIS-up genes - mean of
  ## log2(fold change) across CIS-down genes".'

  .consistencyCheck(nametype, "CISSign")

  sign_df <- CIS_Robertson
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  sign_up <- sign_df[grep("up", sign_df$class), ]
  sign_down <- sign_df[grep("down", sign_df$class), ]

  .percentageOfGenesUsed("CISSign", datasetm, sign_up$SYMBOL, "up")
  .percentageOfGenesUsed("CISSign", datasetm, sign_down$SYMBOL, "down")

  datasetm_n <- datasetm[intersect(row.names(datasetm), sign_df$SYMBOL), ]
  datasetm_n <- log2(datasetm_n + 1) - log2(rowMedians(datasetm_n) + 1)
  score_up <- colMeans(
    datasetm_n[intersect(row.names(datasetm_n), sign_up$SYMBOL), ])
  score_down <- colMeans(
    datasetm_n[intersect(row.names(datasetm_n), sign_down$SYMBOL), ])

  score <- score_up - score_down

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "CIS_Robertson", datasetm))
}


#' Glycolysis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' glycolysisSign(dataset = ovse)
#'
#' @export
glycolysisSign <- function(
    dataset, nametype = "SYMBOL", author = "Zhang",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "glycolysisSign", author)

  sign_df <- get(paste0("Glycolysis_", author))
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- log2(datasetm + 1)

  score <- .coeffScore(sign_df, datasetm_n, "glycolysisSign", author = author)

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("Glycolysis_", author), datasetm))
}


#' Autophagy Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data(ovse)
#' autophagySign(dataset = ovse)
#'
#' @export
autophagySign <- function(
    dataset, nametype = "SYMBOL", author = "Xu",
    whichAssay = "norm_expr", hgReference = "hg38") {

  .consistencyCheck(nametype, "autophagySign", author)

  sign_df <- get(paste0("Autophagy_", author))
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "ChenM") {
    dataset <- .dataTransformation(
      dataset, datasetm, "FPKM", hgReference, nametype)
    datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
    OSscore <- .coeffScore(
      sign_df[sign_df$class == "OS", ], datasetm_n,
      "autophagySign", "OS", author = author)
    DFSscore <- .coeffScore(
      sign_df[sign_df$class == "DFS", ], datasetm_n,
      "autophagySign", "DFS", author = author)
    score <- data.frame(OSscore, DFSscore)
    colnames(score) <- c("Autophagy_ChenM_OS", "Autophagy_ChenM_DFS")
    return(.returnAsInput(
      userdata = dataset, result = t(score), SignName = "", datasetm))
  } else {
    score <- .coeffScore(
      sign_df, datasetm, "autophagySign", author = author)
    return(.returnAsInput(
      userdata = dataset, result = score,
      SignName = paste0("Autophagy_", author), datasetm))
  }
}


#' Extracellular Matrix Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA gsva ssgseaParam
#'
#' @examples
#' data(ovse)
#' ECMSign(dataset = ovse)
#'
#' @export
ECMSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr", ...) {

  .consistencyCheck(nametype, "ECMSign")

  sign_df <- ECM_Chakravarthy
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  sign_up <- sign_df[grep("ECMup", sign_df$class), ]
  sign_down <- sign_df[grep("ECMdown", sign_df$class), ]

  .percentageOfGenesUsed("ECMSign", datasetm, sign_up$SYMBOL, "up")
  .percentageOfGenesUsed("ECMSign", datasetm, sign_down$SYMBOL, "down")

  gene_sets <- list(
    ECM_Chakravarthy_up = sign_up$SYMBOL,
    ECM_Chakravarthy_down = sign_down$SYMBOL)

  dots <- list(...)
  args <- .matchArguments(dots, list(
    exprData = datasetm, geneSets = gene_sets, normalize = FALSE))
  gsvaPar <- do.call(ssgseaParam, args)
  gsva_count <- gsva(gsvaPar, verbose = FALSE)

  return(.returnAsInput(
    userdata = dataset, result = gsva_count, SignName = "", datasetm))
}


#' Homologous Recombination Deficiency Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @importFrom matrixStats rowMedians
#' @importFrom stats t.test
#'
#' @examples
#' data(ovse)
#' HRDSSign(dataset = ovse)
#'
#' @export
HRDSSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "HRDSSign")

  sign_df <- HRDS_Lu
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- datasetm - rowMedians(datasetm)

  .percentageOfGenesUsed("HRDSSign", datasetm, sign_df$SYMBOL)

  HRDS_P <- datasetm_n[intersect(
    row.names(datasetm_n), sign_df[sign_df$coeff == 1, ]$SYMBOL), ]
  HRDS_N <- datasetm_n[intersect(
    row.names(datasetm_n), sign_df[sign_df$coeff == -1, ]$SYMBOL), ]

  score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(x) {
    tmp <- t.test(HRDS_P[, x], HRDS_N[, x], alternative = "two.sided")
    tmp[["statistic"]] }))

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "HRDS_Lu", datasetm))
}


#' Adult Intestinal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' ISCSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
ISCSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "ISCSign")

  sign_df <- ISC_MerlosSuarez
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
  sign_list <- split(sign_df$SYMBOL, sign_df$class)

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- if (inputType == "rnaseq") {
    log2(datasetm + 1)
  } else { datasetm }
  datasetm_n <- datasetm_n - rowMeans(datasetm_n)

  .percentageOfGenesUsed(
    "ISCSign", datasetm_n, sign_list$`ISCEphB2`, "ISCEphB2")
  .percentageOfGenesUsed(
    "ISCSign", datasetm_n, sign_list$`LateTA`, "LateTA")
  .percentageOfGenesUsed(
    "ISCSign", datasetm_n, sign_list$`ISCLgr5`, "ISCLgr5")
  .percentageOfGenesUsed(
    "ISCSign", datasetm_n, sign_list$`Prolif`, "Prolif")

  scores <- vapply(sign_list, function(x) {
    colMeans(datasetm_n[intersect(
      x, row.names(datasetm_n - colMeans(datasetm_n))), ])
  }, double(ncol(datasetm_n)))
  colnames(scores) <- paste0("ISC_MerlosSuarez_", colnames(scores))

  return(.returnAsInput(
    userdata = dataset, result = t(scores), SignName = "", datasetm))
}


#' VEGF Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' VEGFSign(dataset = ovse)
#'
#' @export
VEGFSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "VEGFSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- log2(datasetm + 1)

  score <- .statScore(
    VEGF_Hu$SYMBOL, datasetm_n, nametype, "mean", "VEGFSign")

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "VEGF_Hu", datasetm))
}


#' DNA Repair Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' DNArepSign(dataset = ovse, inputType = "rnaseq")
#'
#' @export
DNArepSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    whichAssay = "norm_expr") {

  .consistencyCheck(nametype, "DNArepSign")

  sign_df <- DNArep_Kang
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm <- .getMatrix(dataset, whichAssay)

  .percentageOfGenesUsed("DNArepSign", datasetm, sign_df$SYMBOL)

  datasetm_n <- if (inputType == "rnaseq") {
    log2(datasetm + 1)
  } else { datasetm }
  datasetm_n <- datasetm_n[row.names(datasetm_n) %in% sign_df$SYMBOL, ]
  datasetm_n <- scale(t(datasetm_n), center = TRUE, scale = FALSE)

  medianexp <- apply(datasetm_n, 2, median)

  genes_h <- intersect(
    colnames(datasetm_n), sign_df[sign_df$class == "high", ]$SYMBOL)
  genes_l <- intersect(
    colnames(datasetm_n), sign_df[sign_df$class == "low", ]$SYMBOL)

  score <- rowSums(datasetm_n[, genes_h] > medianexp[genes_h]) +
    rowSums(datasetm_n[, genes_l] < medianexp[genes_l])

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "DNArep_Kang", datasetm))
}


#' IPSOV Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{ssgseaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA ssgseaParam gsva
#'
#' @examples
#' data(ovse)
#' IPSOVSign(dataset = ovse)
#'
#' @export
IPSOVSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    whichAssay = "norm_expr", ...) {

  .consistencyCheck(nametype, "IPSOVSign")

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- if (inputType == "rnaseq") { log2(datasetm + 1)
  } else { datasetm }
  datasetm_n <- scale(datasetm, center = TRUE, scale = TRUE)

  sign_df <- IPSOV_Shen
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm),]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)

  .percentageOfGenesUsed("IPSOVSign", datasetm, sign_df$SYMBOL)

  dots <- list(...)
  args <- .matchArguments(dots, list(
    exprData = datasetm_n, geneSets = sign_list, normalize = FALSE))
  gsvaPar <- do.call(ssgseaParam, args)
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)

  sign_class <- unique(sign_df[,c(2,3)])
  sign_class <- sign_class[sign_class$class %in% row.names(gsva_matrix), ]
  columnNA <- .managena(datasetm = gsva_matrix, genes = sign_class$class)
  score <- colSums(
    gsva_matrix[sign_class$class, ] * sign_class$coeff, na.rm = TRUE)
  score[columnNA > 0.9] <- NA

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "IPSOV_Shen", datasetm))

  return(dataset)
}


#' EMT-Inflammation Combined Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param weighted logical value, saying whether the score should be calculated
#' with or without weights.
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' CombinedSign(dataset = ovse)
#'
#' @export
CombinedSign <- function(
        dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
        hgReference = "hg38", weighted = FALSE){

  # Combination of EMT_Thompson and Tinflam_Thompson signatures

  .consistencyCheck(nametype, "CombinedSign")

  datasetm <- .getMatrix(dataset, whichAssay)

  res_emt <- EMTSign(
    dataset, nametype = nametype, author = "Thompson", whichAssay = whichAssay,
    hgReference = hgReference)[["EMT_Thompson"]]
  res_inf <- TinflamSign(
    dataset, nametype = nametype, author = "Thompson", whichAssay = whichAssay,
    hgReference = hgReference)[["Tinflam_Thompson"]]

  score <- if (weighted){-0.60*res_emt + 0.19*res_inf} else {res_inf-res_emt}

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Combined_Thompson", datasetm))
}


#' Antigen Processing Machinery Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsvaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA gsva gsvaParam
#'
#' @examples
#' data(ovse)
#' APMSign(dataset = ovse)
#'
#' @export
APMSign <- function(
    dataset, nametype = "SYMBOL", inputType = "microarray",
    author = "Wang", whichAssay = "norm_expr", hgReference = "hg38", ...) {

  .consistencyCheck(nametype, "APMSign", author)

  datasetm <- .getMatrix(dataset, whichAssay)

  if (author == "Wang") {
    sign_df <- APM_Wang
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
    sign_list <- list(sign_df$SYMBOL)
    names(sign_list) <- "APM_Wang"

    .percentageOfGenesUsed("APMSign", datasetm, sign_df$SYMBOL)

    dots <- list(...)
    args <- .matchArguments(dots, list(
      exprData = datasetm, geneSets = sign_list, kcdf = "Gaussian"))

    gsvaPar <- do.call(gsvaParam, args)
    score <- gsva(gsvaPar, verbose = FALSE)
    score <- as.vector(score)

  } else if (author == "Thompson") {

    dataset <- .dataTransformation(
      dataset, datasetm, "CPM", hgReference, nametype)
    datasetm_n <- as.matrix(assays(dataset)[["CPM"]])
    datasetm_n <- log2(datasetm_n+1)

    sign_df <- APM_Thompson
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed("APMSign", datasetm_n, sign_df$SYMBOL)

    z_score <- t(
      scale(t(datasetm_n[intersect(sign_df$SYMBOL, rownames(datasetm_n)), ])))
    score <- colSums(log2(z_score-min(z_score)+1))
  }

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = paste0("APM_", author), datasetm))
}


#' Melanocytic Plasticity Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' MPSSign(dataset = ovse)
#'
#' @export
MPSSign <- function(
        dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
        hgReference = "hg38"){

  .consistencyCheck(nametype, "MPSSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- MPS_PerezGuijarro
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  dataset <- .dataTransformation(
    dataset = dataset, data = datasetm, trans = "FPKM",
    hgReference = hgReference, nametype = nametype)
  datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])

  score <- .coeffScore(
      sdata = sign_df, datasetm = datasetm_n, namesignature = "MPSSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "MPS_PerezGuijarro", datasetm))
}


#' Immune-Related Genes Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' IRGSign(dataset = ovse)
#'
#' @export
IRGSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "IRGSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- IRG_Yang
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  score <- .coeffScore(
    sdata = sign_df, datasetm = datasetm, namesignature = "IRGSign")

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "IRG_Yang", datasetm))
}


#' Pan-Fibroblast TGFB Response Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' TGFBSign(dataset = ovse)
#'
#' @export
TGFBSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "TGFBSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- TGFB_Mariathasan
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed("TGFBSign", datasetm, sign_df$SYMBOL)

  pca<-prcomp(t(datasetm[
    rownames(datasetm) %in% TGFB_Mariathasan$SYMBOL, ]), scale = TRUE)$x
  score <- pca[, "PC1"]

  return(.returnAsInput(
    userdata = dataset, result = score, SignName="TGFB_Mariathasan", datasetm))
}


#' Adenosine Signaling Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsvaParam}}
#' function.
#'
#' @inherit EMTSign return
#'
#' @importFrom GSVA gsvaParam gsva
#'
#' @examples
#' data(ovse)
#' ADOSign(dataset = ovse)
#'
#' @export
ADOSign <- function(
  dataset, nametype = "SYMBOL", whichAssay = "norm_expr", ...){

  .consistencyCheck(nametype, "ADOSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- ADO_Sidders
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)
  sign_list <- list(sign_df$SYMBOL)
  names(sign_list) <- "ADO"

  .percentageOfGenesUsed("ADOSign", datasetm, sign_df$SYMBOL)

  dots <- list(...)
  args <- .matchArguments(dots, list(
    exprData = datasetm, geneSets = sign_list, kcdf = "Poisson"))
  gsvaPar <- do.call(gsvaParam, args)
  score <- as.vector(gsva(gsvaPar, verbose = FALSE))

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "ADO_Sidders", datasetm))
}


#' MITFlow/PTENneg Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' MITFlowPTENnegSign(dataset = ovse)
#'
#' @export
MITFlowPTENnegSign <- function(
  dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "MITFlowPTENnegSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- MITFlowPTENneg_Cabrita
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed("MITFlowPTENnegSign", datasetm, sign_df$SYMBOL)

  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm),]
  datasetm <- log2(datasetm+1)
  score <- cor(
    datasetm[sign_df$SYMBOL,], sign_df$coeff, method = "pearson")

  return(.returnAsInput(
    userdata = dataset, result = as.vector(score),
    SignName = "MITFlowPTENneg_Cabrita", datasetm))
}


#' LRRC15 CAF Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' LRRC15CAFSign(dataset = ovse)
#'
#' @export
LRRC15CAFSign <- function(
  dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "LRRC15CAFSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- LRRC15CAF_Dominguez
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed("LRRC15CAFSign", datasetm, sign_df$SYMBOL)

  datasetm_n <- datasetm[rownames(datasetm) %in% sign_df$SYMBOL,]
  datasetm_n <- log2(datasetm_n + 1)
  pca <- prcomp(t(datasetm_n), center=TRUE, scale=TRUE)
  weights <- pca$rotation[,1]
  contribution <- abs(weights) / sum(abs(weights))
  score <- apply(datasetm_n, 2, function(x) weighted.mean(x, contribution))

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "LRRC15CAF_Dominguez", datasetm))
}



#' Breast Cancer Subtypes Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
SCSubtypeSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "SCSubtypeSign")
  .isMalignantCheck(isMalignant, dataset)

  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm_n <- log2(datasetm + 1)

  sign_df <- SCSubtype_Wu
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "SCSubtypeSign", datasetm_n,
    sign_df$SYMBOL[sign_df$class == "Basal"], "Basal")
  .percentageOfGenesUsed(
    "SCSubtypeSign", datasetm_n,
    sign_df$SYMBOL[sign_df$class == "Her2E"], "Her2E")
  .percentageOfGenesUsed(
    "SCSubtypeSign", datasetm_n,
    sign_df$SYMBOL[sign_df$class == "LumA"], "LumA")
  .percentageOfGenesUsed(
    "SCSubtypeSign", datasetm_n,
    sign_df$SYMBOL[sign_df$class == "LumB"], "LumB")

  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)

  names(sign_list) <- paste0("SCSubtype_Wu_", names(sign_list))

  datasetm_n <- datasetm_n[,isMalignant]

  scores <- as.data.frame(lapply(sign_list, function(x) {
      s <- rep(NA, ncol(datasetm))
      if (length(x)>1) {score <- colMeans(datasetm_n[x,], na.rm = TRUE)
      } else {score <- datasetm_n[x,]}
      s[isMalignant] <- score
      s }))

  return(.returnAsInput(
    userdata = dataset, result = t(scores),
    SignName = "SCSubtype_Wu", datasetm))
}


#' ICB Response Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' ICBResponseSign(dataset = ovse)
#'
#' @export
ICBResponseSign <- function(
        dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "ICBResponseSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- ICBResponse_Chen
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "ICBResponseSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "responder"], "responder")
  .percentageOfGenesUsed(
    "ICBResponseSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "nonresponder"], "nonresponder")

  resp <- sign_df[grep("^responder", sign_df$class),]
  nresp <- sign_df[grep("nonresponder", sign_df$class),]
  score <- list()

  datasetm_r <- datasetm[intersect(rownames(datasetm), resp$SYMBOL), ]
  sup_resp <- colSums(datasetm_r)
  corr_resp <- apply(
      datasetm_r, 1, function(x) cor(x, sup_resp, method = "pearson"))
  score[["ICBResponse_Chen_responder"]] <- colSums(datasetm_r * corr_resp)

  datasetm_n <- datasetm[intersect(rownames(datasetm), nresp$SYMBOL), ]
  sup_nresp <- colSums(datasetm_n)
  corr_nresp <- apply(
      datasetm_n, 1, function(x) cor(x, sup_nresp, method = "pearson"))
  score[["ICBResponse_Chen_nonresponder"]] <- colSums(datasetm_n * corr_nresp)

  return(.returnAsInput(
    userdata = dataset, result = t(as.data.frame(score)),
    SignName = "", datasetm))
}


#' COX-2-associated Inflammatory Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' COXISSign(dataset = ovse)
#'
#' @export
COXISSign <- function(dataset, nametype = "SYMBOL", whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "COXISSign")
  datasetm <- .getMatrix(dataset, whichAssay)
  datasetm <- log2(datasetm + 1)

  sign_df <- COXIS_Bonavita
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed("COXISSign", datasetm, sign_df$SYMBOL)

  pos <- intersect(rownames(datasetm), sign_df[sign_df$class == "CP", "SYMBOL"])
  datasetm_pos <- datasetm[pos,]
  score_pos <- colMeans(datasetm_pos)

  neg <- intersect(rownames(datasetm), sign_df[sign_df$class == "CI", "SYMBOL"])
  datasetm_neg <- datasetm[neg,]
  score_neg <- colMeans(datasetm_neg)

  score <- score_pos / score_neg

  return(.returnAsInput(
    userdata = dataset, result = score, SignName = "COXIS_Bonavita", datasetm))
}


#'Stress Response Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
stressSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "stressSign")
  .isMalignantCheck(isMalignant, dataset)

  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- PanState_Barkley
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "stressSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "Stress"], "Stress")

  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)
  names(sign_list) <- paste0(names(sign_list), "_Barkley")

  datasetm_n <- datasetm[,isMalignant]
  modules = "Stress_Barkley"
  mod <- sign_list[modules]
  score <- .barkleyFun(
      dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Stress_Barkley", datasetm))
}

#'Interferon Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
interferonSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "interferonSign")

  datasetm <- .getMatrix(dataset, whichAssay)

  .isMalignantCheck(isMalignant, dataset)

  sign_df <- PanState_Barkley
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "interferonSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "Interferon"], "Interferon")

  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)
  names(sign_list) <- paste0(names(sign_list), "_Barkley")

  datasetm_n <- datasetm[,isMalignant]
  modules = "Interferon_Barkley"
  mod <- sign_list[modules]
  score <- .barkleyFun(
      dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Interferon_Barkley", datasetm))
}

#'Oxidative Phosphorylation Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
oxphosSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "oxphosSign")

  datasetm <- .getMatrix(dataset, whichAssay)

  .isMalignantCheck(isMalignant, dataset)

  sign_df <- PanState_Barkley
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "oxphosSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "Oxphos"], "Oxphos")


  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)
  names(sign_list) <- paste0(names(sign_list), "_Barkley")

  datasetm_n <- datasetm[,isMalignant]
  modules = "Oxphos_Barkley"
  mod <- sign_list[modules]
  score <- .barkleyFun(
      dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Oxphos_Barkley", datasetm))
}

#'Metal Response Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
metalSign <- function(
    dataset, nametype = "SYMBOL", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "metalSign")

  datasetm <- .getMatrix(dataset, whichAssay)

  .isMalignantCheck(isMalignant, dataset)

  sign_df <- PanState_Barkley
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  .percentageOfGenesUsed(
    "metalSign", datasetm,
    sign_df$SYMBOL[sign_df$class == "Metal"], "Metal")


  sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
  sign_list <- split(sign_df$SYMBOL, sign_df$class)
  names(sign_list) <- paste0(names(sign_list), "_Barkley")

  datasetm_n <- datasetm[,isMalignant]
  modules = "Metal_Barkley"
  mod <- sign_list[modules]
  score <- .barkleyFun(
      dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "Metal_Barkley", datasetm))
}

#' Cellular States Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#'
#' @export
stateSign <- function(
    dataset, nametype = "SYMBOL", author = "Barkley", whichAssay = "norm_expr",
    isMalignant = NULL, hgReference = "hg38") {

  .consistencyCheck(nametype, "stateSign")
  .isMalignantCheck(isMalignant, dataset)

  if (author == "Neftel") {

    if(nrow(dataset)<3000){stop(
      "dataset must have at least 3000 genes to compute the signature")}

    datasetm <- .getMatrix(dataset, whichAssay)
    dataset <- .dataTransformation(
      dataset, datasetm, "TPM", hgReference, nametype)
    datasetm_n <- as.matrix(assays(dataset)[["TPM"]])

    sign_df <- State_Neftel
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "MES2"], "MES2")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "MES1"], "MES1")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "AC"], "AC")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "OPC"], "OPC")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "NPC1"], "NPC1")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "NPC2"], "NPC2")

    sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
    sign_list <- split(sign_df$SYMBOL, sign_df$class)
    names(sign_list) <- paste0("State_Neftel_", names(sign_list))

    datasetm_n <- datasetm_n[,isMalignant]
    exp_lev <- log2(datasetm_n/10+1)
    rel_exp <- exp_lev - rowMeans(exp_lev,  na.rm = TRUE)

    agg_exp <- log2(rowMeans(datasetm_n, na.rm = TRUE)+1)
    ea_bin <- split(
      sort(agg_exp, na.last = TRUE), factor(
        sort(round(x = rank(agg_exp) %% 30, digits = 0))))
    ea_bin <- lapply(ea_bin, function(x){names(x)})

    scores <- as.data.frame(lapply(sign_list, function(x){
      Gcont <- unlist(lapply(x, function(y){
        u <- NULL
        for (i in seq_along(ea_bin)) {
          if (y %in% ea_bin[[i]]) {
            u <- i
            break}}
        sample(ea_bin[[u]][!(ea_bin[[u]] %in% x)], 100)}))
      score <- rep(NA, ncol(dataset))
      SC <- colMeans(
        rel_exp[x,], na.rm = TRUE)-colMeans(rel_exp[Gcont,], na.rm = TRUE)
      score[isMalignant] <- SC
      score
    }))

    return(.returnAsInput(
      userdata = dataset, result = t(scores), SignName = "", datasetm))

  } else if (author == "Barkley") {

    datasetm <- .getMatrix(dataset, whichAssay)
    sign_df <- PanState_Barkley
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Alveolar"], "Alveolar")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Basal"], "Basal")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Squamous"], "Squamous")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Glandular"], "Glandular")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "Ciliated"], "Ciliated")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "AC"], "AC")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "OPC"], "OPC")
    .percentageOfGenesUsed(
      "stateSign", datasetm,
      sign_df$SYMBOL[sign_df$class == "NPC"], "NPC")


    sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm), ]
    sign_list <- split(sign_df$SYMBOL, sign_df$class)
    names(sign_list) <- paste0("State_Barkley_", names(sign_list))

    datasetm_n <- datasetm[,isMalignant]
    modules = c(
        "State_Barkley_Alveolar", "State_Barkley_Basal",
        "State_Barkley_Squamous", "State_Barkley_Glandular",
        "State_Barkley_Ciliated", "State_Barkley_AC",
        "State_Barkley_OPC", "State_Barkley_NPC")
    mod <- sign_list[modules]
    score <- .barkleyFun(
        dataset = datasetm_n, mod = mod, n = ncol(datasetm), isMalignant)

    return(.returnAsInput(
      userdata = dataset, result = score, SignName = "", datasetm))
  } else if (author == "Tirosh") {

    if(nrow(dataset)<2500){stop(
      "dataset must have at least 2500 genes to compute the signature")}

    datasetm <- .getMatrix(dataset, whichAssay)
    dataset <- .dataTransformation(
      dataset, datasetm, "TPM", hgReference, nametype)
    datasetm_n <- as.matrix(assays(dataset)[["TPM"]])

    sign_df <- State_Tirosh
    sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "MITF"], "MITF")
    .percentageOfGenesUsed(
      "stateSign", datasetm_n,
      sign_df$SYMBOL[sign_df$class == "AXL"], "AXL")

    sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
    sign_list <- split(sign_df$SYMBOL, sign_df$class)
    names(sign_list) <- paste0("State_Tirosh_", names(sign_list))

    datasetm_n <- datasetm_n[,isMalignant]
    exp_lev <- log2(datasetm_n/10+1)
    rel_exp <- exp_lev - rowMeans(exp_lev,  na.rm = TRUE)

    agg_exp <- log2(rowMeans(datasetm_n, na.rm = TRUE)+1)
    ea_bin <- split(
      sort(agg_exp, na.last = TRUE), factor(
        sort(round(x = rank(agg_exp) %% 25, digits = 0))))
    ea_bin <- lapply(ea_bin, function(x){names(x)})

    scores <- as.data.frame(lapply(sign_list, function(x){
      Gcont <- unlist(lapply(x, function(y){
        u <- NULL
        for (i in seq_along(ea_bin)) {
          if (y %in% ea_bin[[i]]) {
            u <- i
            break}}
        sample(ea_bin[[u]][!(ea_bin[[u]] %in% x)], 100)}))
      score <- rep(NA, ncol(dataset))
      SC <- colMeans(
        rel_exp[x,], na.rm = TRUE)-colMeans(rel_exp[Gcont,], na.rm = TRUE)
      score[isMalignant] <- SC
      score
    }))

    return(.returnAsInput(
      userdata = dataset, result = t(scores), SignName = "", datasetm))
  }

}


#' CD39+ CD8+ T cells Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @inherit EMTSign return
#'
#' @examples
#' data(ovse)
#' CD39CD8TcellSign(dataset = ovse)
#'
#' @export
CD39CD8TcellSign <- function(dataset, nametype = "SYMBOL",
                             whichAssay = "norm_expr"){

  .consistencyCheck(nametype, "CD39CD8TcellSign")
  datasetm <- .getMatrix(dataset, whichAssay)

  sign_df <- CD39CD8Tcell_Chow
  sign_df$SYMBOL <- .geneIDtrans(nametype, sign_df$SYMBOL)

  datasetm_n <- t(scale(t(datasetm)))

  score <- .statScore(genes = sign_df$SYMBOL, datasetm = datasetm_n,
                      nametype = nametype, typeofstat = "median",
                      namesignature = "CD39CD8TcellSign")

  return(.returnAsInput(
    userdata = dataset, result = score,
    SignName = "CD39CD8Tcell_Chow", datasetm))
}
