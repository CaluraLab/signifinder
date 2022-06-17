

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
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}}
#' function.
#'
#' @return A SummarizedExperiment object in which the Epithelial and Mesenchymal
#' scores are added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom stats prcomp
#'
#' @examples
#' EMTSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
EMTSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "microarray",
             author = "Miow",
             pvalues = FALSE,
             nperm = 100,
             ...) {
        consistencyCheck(nametype, "EMTSign", author)

        datasetm <- getMatrix(dataset)

        if (author == "Miow") {
            sign_df <- EMT_Miow
            sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

            EL <- sign_df[grep("Epithelial-like", sign_df$class), ]
            ML <- sign_df[grep("Mesenchymal-like", sign_df$class), ]

            percentageOfGenesUsed(
                "EMTSign", datasetm, EL$SYMBOL, "epithelial",
                author = author
            )
            percentageOfGenesUsed(
                "EMTSign", datasetm, ML$SYMBOL, "mesenchymal",
                author = author
            )

            gene_sets <-
                list(
                    Epithelial = EL$SYMBOL,
                    Mesenchymal = ML$SYMBOL
                )
            names(gene_sets) <- paste0("EMT_Miow_", names(gene_sets))

            dots <- list(...)
            kcdftype <-
                ifelse(inputType == "microarray", "Gaussian", "Poisson")
            args <- matchArguments(
                dots,
                list(
                    expr = datasetm,
                    gset.idx.list = gene_sets,
                    method = "ssgsea",
                    kcdf = kcdftype,
                    ssgsea.norm = FALSE,
                    verbose = FALSE
                )
            )
            gsva_matrix <- do.call(gsva, args)

            if (pvalues) {
                gsva_pval <- GSVAPvalues(
                    expr = datasetm,
                    gset.idx.list = gene_sets,
                    gsvaResult = gsva_matrix,
                    nperm = nperm,
                    args = args
                )
                gsva_matrix <- rbind(gsva_matrix, gsva_pval)
            }

            return(
                returnAsInput(
                    userdata = dataset,
                    result = gsva_matrix,
                    SignName = "",
                    datasetm
                )
            )
        } else {
            if (author == "Mak") {
                sign_df <- EMT_Mak
                sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

                Sign_E <- sign_df$SYMBOL[sign_df$class == "E"]
                Sign_M <- sign_df$SYMBOL[sign_df$class == "M"]

                percentageOfGenesUsed(
                    "EMTSign", datasetm, Sign_E, "epithelial",
                    author = author
                )
                percentageOfGenesUsed(
                    "EMTSign", datasetm, Sign_M, "mesenchymal",
                    author = author
                )

                columnNA <- managena(datasetm, genes = sign_df$SYMBOL)
                score <- colMeans(
                    datasetm[intersect(Sign_M, row.names(datasetm)), ]
                ) -
                    colMeans(datasetm[intersect(Sign_E, row.names(datasetm)), ])
                score[columnNA > 0.9] <- NA
            } else if (author == "Cheng") {
                sign_df <- EMT_Cheng
                sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

                percentageOfGenesUsed(
                    "EMTSign", datasetm, sign_df$SYMBOL,
                    author = author
                )

                datasetm_n <- datasetm[
                    intersect(row.names(datasetm), sign_df$SYMBOL),
                ]
                datasetm_n <-
                    if (inputType == "rnaseq") {
                        log2(datasetm_n + 0.01)
                    } else {
                        datasetm_n
                    }
                columnNA <- managena(
                    datasetm = datasetm_n,
                    genes = sign_df$SYMBOL
                )
                score <- prcomp(t(datasetm_n))$x[, 1]
                score[columnNA > 0.9] <- NA
            }
            return(returnAsInput(
                userdata = dataset,
                result = score,
                SignName = paste0("EMT_", author),
                datasetm
            ))
        }
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
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' pyroptosisSign(dataset = subovse)
#'
#' @export
pyroptosisSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "rnaseq",
             author = "Ye",
             hgReference = "hg38") {
        consistencyCheck(nametype, "pyroptosisSign", author)

        sign_df <- get(paste0("Pyroptosis_", author))
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)

        if (author == "Ye") {
            dataset <- dataTransformation(
                dataset, datasetm, "FPKM", hgReference, nametype
            )
            datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
            datasetm_n <- scale(datasetm_n)
        } else if (author == "Shao") {
            if (inputType == "rnaseq") {
                dataset <- dataTransformation(
                    dataset, datasetm, "FPKM", hgReference, nametype
                )
                datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
            } else {
                datasetm_n <- datasetm
            }
        } else if (author == "Lin") {
            dataset <- dataTransformation(
                dataset, datasetm, "TPM", hgReference, nametype
            )
            datasetm_n <- as.matrix(assays(dataset)[["TPM"]])
        } else {
            datasetm_n <- datasetm
        }

        score <-
            coeffScore(sign_df, datasetm_n, "pyroptosisSign", author = author)

        return(returnAsInput(
            userdata = dataset,
            result = score,
            SignName = paste0("Pyroptosis_", author),
            datasetm
        ))
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
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' ferroptosisSign(dataset = subovse)
#'
#' @export
ferroptosisSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "rnaseq",
             author = "Ye",
             hgReference = "hg38") {
        consistencyCheck(nametype, "ferroptosisSign", author)

        sign_df <- get(paste0("Ferroptosis_", author))
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)

        if (author %in% c("Liu", "Ye") & inputType == "rnaseq") {
            dataset <- dataTransformation(
                dataset, datasetm, "FPKM", hgReference, nametype
            )
            datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
            datasetm_n <-
                datasetm_n[rownames(datasetm_n) %in% sign_df$SYMBOL, ]
        } else {
            datasetm_n <- datasetm[rownames(datasetm) %in% sign_df$SYMBOL, ]
        }

        score <-
            coeffScore(sign_df, datasetm_n, "ferroptosisSign", author = author)

        if (author == "Liang") {
            score <- exp(score)
        }

        return(returnAsInput(
            userdata = dataset,
            result = score,
            SignName = paste0("Ferroptosis_", author),
            datasetm
        ))
    }


#' Lipid Metabolism Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Lipid scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' lipidMetabolismSign(dataset = subovse)
#'
#' @export
lipidMetabolismSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "lipidMetabolismSign")

    sign_df <- LipidMetabolism_Zheng
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)
    score <- coeffScore(sign_df, datasetm, "lipidMetabolismSign")

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "LipidMetabolism_Zheng",
            datasetm
        )
    )
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
#' @examples
#' hypoxiaSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
hypoxiaSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "microarray") {
        consistencyCheck(nametype, "hypoxiaSign")

        datasetm <- getMatrix(dataset)
        datasetm_n <-
            if (inputType == "rnaseq") {
                log2(datasetm + 0.01)
            } else {
                datasetm
            }
        score <-
            statScore(
                Hypoxia_Buffa$SYMBOL,
                datasetm = abs(datasetm_n),
                nametype = nametype,
                typeofstat = "median",
                namesignature = "hypoxiaSign"
            )

        return(
            returnAsInput(
                userdata = dataset,
                result = as.vector(scale(score)),
                SignName = "Hypoxia_Buffa",
                datasetm
            )
        )
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
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' immunoScoreSign(dataset = subovse)
#'
#' @export
immunoScoreSign <-
    function(dataset,
             nametype = "SYMBOL",
             author = "Hao",
             inputType = "rnaseq",
             hgReference = "hg38") {
        consistencyCheck(nametype, "immunoScoreSign", author)

        datasetm <- getMatrix(dataset)

        if (author == "Hao") {
            sign_df <- ImmunoScore_Hao
            sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

            g <- intersect(row.names(datasetm), sign_df$SYMBOL)

            percentageOfGenesUsed(
                "immunoScoreSign", datasetm, sign_df$SYMBOL,
                author = author
            )

            if (inputType == "rnaseq") {
                dataset <- dataTransformation(
                    dataset, datasetm, "FPKM", hgReference, nametype
                )
                datasetm_n <-
                    log2(as.matrix(assays(dataset)[["FPKM"]]) + 0.01)
                datasetm_n <- datasetm_n[g, ]
            } else {
                datasetm_n <- datasetm[g, ]
            }

            sign_df <- sign_df[sign_df$SYMBOL %in% g, ]
            columnNA <- managena(datasetm_n, g)
            SE <- (sign_df$HR - sign_df$`95CI_L`) / 1.96
            k <- (1 - sign_df$HR) / SE
            score <-
                unlist(lapply(seq_len(ncol(datasetm_n)), function(p) {
                    sum(k * datasetm_n[, p], na.rm = TRUE)
                }))
            score[columnNA > 0.9] <- NA
        } else if (author == "Roh") {
            datasetm_n <- log2(datasetm + 1)
            score <- statScore(
                ImmunoScore_Roh$SYMBOL,
                datasetm = datasetm_n,
                nametype = nametype,
                typeofstat = "meang",
                namesignature = "immunoScoreSign",
                author = author
            )
        }

        return(returnAsInput(
            userdata = dataset,
            result = score,
            SignName = paste0("ImmunoScore_", author),
            datasetm
        ))
    }


#' ConsensusOV Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
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
#' @examples
#' consensusOVSign(dataset = subovse)
#'
#' @export
consensusOVSign <- function(dataset,
                            nametype = "SYMBOL",
                            method = "consensusOV",
                            ...) {
    consistencyCheck(nametype, "consensusOVSign")

    datasetm <- getMatrix(dataset)

    if (nametype != "ENTREZID") {
        genename <- mapIds(
            org.Hs.eg.db,
            keys = row.names(datasetm),
            column = "ENTREZID",
            keytype = nametype,
            multiVals = "first"
        )
        datasetm_n <- datasetm[!is.na(genename), ]
        genename <- genename[!is.na(genename)]
        datasetm_n <- datasetm_n[!duplicated(genename), ]
        genename <- genename[!duplicated(genename)]
    } else {
        genename <- row.names(datasetm)
        datasetm_n <- datasetm
    }

    consensus_subtypes <-
        get.subtypes(
            expression.dataset = datasetm_n,
            entrez.ids = genename,
            method = method,
            ...
        )
    scores <- consensus_subtypes$rf.probs
    colnames(scores) <- paste0(
        "ConsensusOV_Chen_", substring(colnames(scores), 1, 3)
    )

    return(returnAsInput(
        userdata = dataset,
        result = t(scores),
        SignName = "",
        datasetm
    ))
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
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' IPSSign(dataset = subovse)
#'
#' @export
IPSSign <-
    function(dataset,
             nametype = "SYMBOL",
             hgReference = "hg38") {
        ## code adapted from https://github.com/icbi-lab/Immunophenogram

        consistencyCheck(nametype, "IPSSign")

        sign_df <- IPS_Charoentong
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)
        sample_names <- colnames(datasetm)

        percentageOfGenesUsed("IPSSign", datasetm, sign_df$SYMBOL)

        MISSING_GENES <- sign_df$SYMBOL[
            is.na(match(sign_df$SYMBOL, rownames(datasetm)))
        ]
        if (length(MISSING_GENES) > 0) {
            message(
                "Differently named or missing genes: ",
                MISSING_GENES,
                "\n"
            )
        }

        dataset <- dataTransformation(
            dataset, datasetm, "TPM", hgReference, nametype
        )
        datasetm_n <- log2(as.matrix(assays(dataset)[["TPM"]]) + 1)
        datasetm_n <-
            datasetm_n[rownames(datasetm_n) %in% sign_df$SYMBOL, ]

        IPS <-
            NULL
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
                WEIGHT[k] <-
                    mean(sign_df$WEIGHT[which(sign_df$NAME == gen)])
                k <- k + 1
            }
            WG <- MIG * WEIGHT
            MHC[i] <- mean(WG[which(
                unique(sign_df$NAME) %in% unique(
                    sign_df$NAME[sign_df$class == "MHC"]
                )
            )], na.rm = TRUE)
            CP[i] <- mean(WG[which(
                unique(sign_df$NAME) %in% unique(
                    sign_df$NAME[sign_df$class == "CP"]
                )
            )], na.rm = TRUE)
            EC[i] <- mean(WG[which(
                unique(sign_df$NAME) %in% unique(
                    sign_df$NAME[sign_df$class == "EC"]
                )
            )], na.rm = TRUE)
            SC[i] <- mean(WG[which(
                unique(sign_df$NAME) %in% unique(
                    sign_df$NAME[sign_df$class == "SC"]
                )
            )], na.rm = TRUE)
            AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
            IPS[i] <- ipsmap(AZ[i])
        }

        ipsres <- data.frame(IPS, MHC, CP, EC, SC)
        rownames(ipsres) <- sample_names
        colnames(ipsres) <- c(
            "IPS_Charoentong",
            "IPS_Charoentong_MHC",
            "IPS_Charoentong_CP",
            "IPS_Charoentong_EC",
            "IPS_Charoentong_SC"
        )
        return(returnAsInput(
            userdata = dataset,
            result = t(ipsres),
            SignName = "",
            datasetm
        ))
    }


#' Core Matrisome Gene signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' matrisomeSign(dataset = subovse)
#'
#' @export
matrisomeSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "matrisomeSign")

    datasetm <- getMatrix(dataset)
    score <-
        statScore(
            Matrisome_Yuzhalin$SYMBOL,
            datasetm,
            nametype,
            "median",
            "matrisomeSign"
        )

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "Matrisome_Yuzhalin",
            datasetm
        )
    )
}


#' Mitotic Index
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' mitoticIndexSign(dataset = subovse)
#'
#' @export
mitoticIndexSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "mitoticIndexSign")

    datasetm <- getMatrix(dataset)

    score <- statScore(
        MitoticIndex_Yang$SYMBOL,
        datasetm,
        nametype,
        "mean",
        "mitoticIndexSign"
    )

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "MitoticIndex_Yang",
            datasetm
        )
    )
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
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' immuneCytSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
immuneCytSign <- function(dataset,
                          nametype = "SYMBOL",
                          inputType = "microarray",
                          author = "Rooney",
                          hgReference = "hg38") {
    consistencyCheck(nametype, "immuneCytSign", author)

    datasetm <- getMatrix(dataset)
    if (author == "Rooney") {
        if (inputType == "rnaseq") {
            dataset <- dataTransformation(
                dataset, datasetm, "TPM", hgReference, nametype
            )
            datasetm_n <- as.matrix(assays(dataset)[["TPM"]]) + 0.01
        } else {
            datasetm_n <- datasetm
        }
        score <- statScore(
            ImmuneCyt_Rooney$SYMBOL,
            datasetm_n,
            nametype,
            "meang",
            "ImmuneCytSign",
            author = author
        )
    } else if (author == "Davoli") {
        datasetm_n <-
            datasetm[row.names(datasetm) %in% ImmuneCyt_Davoli$SYMBOL, ]
        datasetm_r <- apply(datasetm_n, 1, rank)
        datasetm_r <- (datasetm_r - 1) / (nrow(datasetm_r) - 1)
        score <- statScore(
            ImmuneCyt_Davoli$SYMBOL,
            t(datasetm_r),
            nametype,
            "meang",
            "ImmuneCytSign",
            author = author
        )
    }
    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = paste0("ImmuneCyt_", author),
        datasetm
    ))
}


#' IFN-gamma Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' IFNSign(dataset = subovse)
#'
#' @export
IFNSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "IFNSign")

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm + 0.01)
    score <- statScore(
        IFN_Ayers$SYMBOL, datasetm_n,
        nametype, "mean", "IFNSign"
    )

    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = "IFN_Ayers",
        datasetm
    ))
}


#' ExpandedImmune Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' expandedImmuneSign(dataset = subovse)
#'
#' @export
expandedImmuneSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "expandedImmuneSign")

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm + 0.01)
    score <- statScore(
        ExpandedImmune_Ayers$SYMBOL,
        datasetm_n,
        nametype,
        "mean",
        "expandedImmuneSign"
    )

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "ExpandedImmune_Ayers",
            datasetm
        )
    )
}


#' TinflamSign Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' TinflamSign(dataset = subovse)
#'
#' @export
TinflamSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "TinflamSign")

    sign_df <- Tinflam_Ayers
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm + 0.01)

    housekeeping <- intersect(
        row.names(datasetm_n), sign_df$SYMBOL[sign_df$class == "Housekeeping"]
    )
    genes_pred <- intersect(
        row.names(datasetm_n), sign_df$SYMBOL[sign_df$class == "TInflam"]
    )

    housekeeping_m <- apply(datasetm_n[housekeeping, ], 2, mean)
    datasetm_n <-
        sweep(datasetm_n[genes_pred, ], 2, housekeeping_m, FUN = "-")
    score <- coeffScore(
        sign_df[sign_df$SYMBOL %in% genes_pred, ], datasetm_n, "TinflamSign"
    )

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "Tinflam_Ayers",
            datasetm
        )
    )
}


#' Tertiary Lymphoid Structures (TLS) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' TLSSign(dataset = subovse)
#'
#' @export
TLSSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "TLSSign")

    datasetm <- getMatrix(dataset)
    score <- statScore(
        TLS_Cabrita$SYMBOL, datasetm,
        nametype, "mean", "TLSSign"
    )

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "TLS_Cabrita",
            datasetm
        )
    )
}


#' CD49fHi Basal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' stemCellCD49fSign(dataset = subovse)
#'
#' @export
stemCellCD49fSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "StemCellCD49fSign")

    sign_df <- StemCellCD49f_Smith
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm + 0.01)

    score <- coeffScore(sign_df, datasetm_n, "StemCellCD49fSign")

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "StemCellCD49f_Smith",
            datasetm
        )
    )
}


#' Chromosomal instability Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' CINSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
CINSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "microarray") {
        consistencyCheck(nametype, "CINSign")

        datasetm <- getMatrix(dataset)
        datasetm_n <-
            if (inputType == "rnaseq") {
                log2(datasetm + 0.01)
            } else {
                datasetm
            }
        datasetm_n <- scale(datasetm_n, center = TRUE, scale = FALSE)

        score25 <-
            statScore(
                CIN_Carter$SYMBOL[CIN_Carter$class == "CIN25"],
                datasetm_n, nametype, "sum", "CINSign"
            )
        score70 <- statScore(
            CIN_Carter$SYMBOL, datasetm_n,
            nametype, "sum", "CINSign"
        )
        score <- data.frame(score25, score70)
        colnames(score) <- c("CIN_Carter_25", "CIN_Carter_70")

        return(returnAsInput(
            userdata = dataset,
            result = t(score),
            SignName = "",
            datasetm
        ))
    }


#' Cell-cycle Signature classifier
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' cellCycleSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
cellCycleSign <- function(dataset,
                          nametype = "SYMBOL",
                          author = "Lundberg",
                          inputType = "microarray") {
    consistencyCheck(nametype, "cellCycleSign", author)

    datasetm <- getMatrix(dataset)

    if (author == "Lundberg") {
        datasetm_n <- if (inputType == "rnaseq") {
            log2(datasetm + 0.01)
        } else {
            datasetm
        }
        score <- statScore(
            CellCycle_Lundberg$SYMBOL,
            datasetm_n,
            nametype,
            "sum",
            "cellCycleSign",
            author = author
        )
    } else if (author == "Davoli") {
        datasetm_n <-
            datasetm[row.names(datasetm) %in% CellCycle_Davoli$SYMBOL, ]
        datasetm_r <- apply(datasetm_n, 1, rank)
        datasetm_r <- (datasetm_r - 1) / (nrow(datasetm_r) - 1)
        score <- statScore(
            CellCycle_Davoli$SYMBOL,
            t(datasetm_r),
            nametype,
            "meang",
            "cellCycleSign",
            author = author
        )
    }

    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = paste0("CellCycle_", author),
        datasetm
    ))
}


#' Chemokine Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom stats prcomp
#'
#' @examples
#' chemokineSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
chemokineSign <- function(dataset,
                          nametype = "SYMBOL",
                          inputType = "microarray") {
    consistencyCheck(nametype, "chemokineSign")

    sign_df <- Chemokines_Messina
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ChemokineSign", datasetm, sign_df$SYMBOL)

    datasetm_n <-
        datasetm[intersect(row.names(datasetm), sign_df$SYMBOL), ]
    datasetm_n <- if (inputType == "rnaseq") {
        log2(datasetm_n + 0.01)
    } else {
        datasetm_n
    }
    columnNA <- managena(datasetm_n, sign_df$SYMBOL)
    score <-
        prcomp(t(datasetm_n), center = TRUE, scale = TRUE)$x[, 1]
    score[columnNA > 0.9] <- NA

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "Chemokines_Messina",
            datasetm
        )
    )
}

#' Adult Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' ASCSign(dataset = subovse)
#'
#' @export
ASCSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "ASCSign")

    sign_df <- ASC_Smith
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)

    percentageOfGenesUsed("ASCSign", datasetm, sign_df$SYMBOL)

    datasetm_n <-
        log2(datasetm[row.names(datasetm) %in% sign_df$SYMBOL, ] + 1)

    columnNA <- managena(datasetm_n, sign_df$SYMBOL)

    score <- rowSums(
        scale(t(datasetm_n), center = TRUE, scale = TRUE),
        na.rm = TRUE
    )
    score[columnNA > 0.9] <- NA

    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = "ASC_Smith",
        datasetm
    ))
}


#' passON Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param hgReference Human reference genome: "hg19" or "hg38"
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}}
#' function.
#'
#' @return A SummarizedExperiment object in which the passON score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom SummarizedExperiment assays
#' @importFrom stats weighted.mean
#'
#' @examples
#' PassONSign(dataset = subovse)
#'
#' @export
PassONSign <-
    function(dataset,
             nametype = "SYMBOL",
             hgReference = "hg38",
             ...) {
        consistencyCheck(nametype, "PassONSign")

        datasetm <- getMatrix(dataset)
        dataset <- dataTransformation(
            dataset, datasetm, "TPM", hgReference, nametype
        )
        datasetm_n <- as.matrix(assays(dataset)[["TPM"]])

        sign_df <- PassON_Du
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)
        sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
        sign_list <- split(sign_df$SYMBOL, sign_df$class)

        percentageOfGenesUsed("PassONSign", datasetm_n, sign_df$SYMBOL)

        dots <- list(...)
        args <- matchArguments(
            dots,
            list(
                expr = datasetm_n,
                gset.idx.list = sign_list,
                method = "ssgsea",
                kcdf = "Poisson",
                ssgsea.norm = TRUE,
                verbose = FALSE
            )
        )

        gsva_matrix <- do.call(gsva, args)

        gsva_mean <- sapply(seq_len(ncol(gsva_matrix)), function(x) {
            weighted.mean(gsva_matrix[, x], lengths(sign_list))
        })

        return(
            returnAsInput(
                userdata = dataset,
                result = gsva_mean,
                SignName = "PassON_Du",
                datasetm
            )
        )
    }

#' IPRES Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param hgReference Human reference genome: "hg19" or "hg38"
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}}
#' function.
#'
#' @return A SummarizedExperiment object in which the passON score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' IPRESSign(dataset = subovse)
#'
#' @export
IPRESSign <-
    function(dataset,
             nametype = "SYMBOL",
             hgReference = "hg38",
             ...) {
        consistencyCheck(nametype, "IPRESSign")

        datasetm <- getMatrix(dataset)
        dataset <- dataTransformation(
            dataset, datasetm, "CPM", hgReference, nametype
        )
        datasetm_n <- log2(as.matrix(assays(dataset)[["CPM"]]))

        sign_df <- IPRES_Hugo
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)
        sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm_n), ]
        sign_list <- split(sign_df$SYMBOL, sign_df$class)

        percentageOfGenesUsed("IPRESSign", datasetm_n, sign_df$SYMBOL)

        dots <- list(...)
        args <- matchArguments(
            dots,
            list(
                expr = datasetm_n,
                gset.idx.list = sign_list,
                method = "ssgsea",
                kcdf = "Gaussian",
                ssgsea.norm = TRUE,
                verbose = FALSE
            )
        )

        gsva_matrix <- do.call(gsva, args)

        score <- rowMeans(sapply(as.data.frame(t(gsva_matrix)), scale))

        return(
            returnAsInput(
                userdata = dataset,
                result = score,
                SignName = "IPRES_Hugo",
                datasetm
            )
        )
    }

#'  CIS (carcinoma-in situ) Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the CIS score
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' CISSign(dataset = subovse)
#'
#' @export
CISSign <- function(dataset, nametype = "SYMBOL") {
    ## code instruction from an e-mail by Jaegil Kim:
    ## 'We first compute "log2(fold change) = log2(expression) - log2(median
    ## expression across samples)" for each gene and define the signature
    ## score as "mean of log2(fold change) across CIS-up genes - mean of
    ## log2(fold change) across CIS-down genes".'

    consistencyCheck(nametype, "CISSign")

    sign_df <- CIS_Robertson
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)

    sign_up <- sign_df[grep("up", sign_df$class), ]
    sign_down <- sign_df[grep("down", sign_df$class), ]

    percentageOfGenesUsed("CISSign", datasetm, sign_up$SYMBOL, "up")
    percentageOfGenesUsed("CISSign", datasetm, sign_down$SYMBOL, "down")

    datasetm_n <-
        datasetm[intersect(row.names(datasetm), sign_df$SYMBOL), ]
    datasetm_n <-
        log2(datasetm_n + 0.01) - log2(rowMedians(datasetm_n) + 0.01)
    score_up <- colMeans(
        datasetm_n[intersect(row.names(datasetm_n), sign_up$SYMBOL), ])
    score_down <- colMeans(
        datasetm_n[intersect(row.names(datasetm_n), sign_down$SYMBOL), ])

    score <- score_up - score_down

    return(
        returnAsInput(
            userdata = dataset,
            result = score,
            SignName = "CIS_Robertson",
            datasetm
        )
    )
}

#' Glycolysis Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the glycolysis score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' glycolysisSign(dataset = subovse)
#'
#' @export
glycolysisSign <-
    function(dataset,
             nametype = "SYMBOL",
             author = "Zhang") {
        consistencyCheck(nametype, "glycolysisSign", author)

        sign_df <- get(paste0("Glycolysis_", author))
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)
        datasetm_n <- log2(datasetm + 0.01)

        score <-
            coeffScore(sign_df, datasetm_n, "glycolysisSign", author = author)

        return(returnAsInput(
            userdata = dataset,
            result = score,
            SignName = paste0("Glycolysis_", author),
            datasetm
        ))
    }

#' Autophagy Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param hgReference Human reference genome: "hg19" or "hg38"
#'
#' @return A SummarizedExperiment object in which the Autophagy score
#' is added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' autophagySign(dataset = subovse)
#'
#' @export
autophagySign <-
    function(dataset,
             nametype = "SYMBOL",
             author = "Xu",
             hgReference = "hg38") {
        consistencyCheck(nametype, "autophagySign", author)

        sign_df <- get(paste0("Autophagy_", author))
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)

        if (author == "ChenM") {
            dataset <- dataTransformation(
                dataset, datasetm, "FPKM", hgReference, nametype
            )
            datasetm_n <- as.matrix(assays(dataset)[["FPKM"]])
            OSscore <- coeffScore(sign_df[sign_df$class == "OS", ], datasetm_n,
                "autophagySign", "OS",
                author = author
            )
            DFSscore <- coeffScore(sign_df[sign_df$class == "DFS", ], datasetm_n,
                "autophagySign", "DFS",
                author = author
            )
            score <- data.frame(OSscore, DFSscore)
            colnames(score) <-
                c("Autophagy_ChenM_OS", "Autophagy_ChenM_DFS")
            return(returnAsInput(
                userdata = dataset,
                result = t(score),
                SignName = "",
                datasetm
            ))
        } else {
            score <-
                coeffScore(sign_df, datasetm, "autophagySign", author = author)
            return(returnAsInput(
                userdata = dataset,
                result = score,
                SignName = paste0("Autophagy_", author),
                datasetm
            ))
        }
    }

#' Extracellular Matrix Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#' @param ... other arguments passed on to the \code{\link[GSVA]{gsva}}
#' function.
#'
#' @return A SummarizedExperiment object in which the ECM scores
#' will be added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom GSVA gsva
#'
#' @examples
#' ECMSign(dataset = subovse)
#'
#' @export
ECMSign <- function(dataset,
                    nametype = "SYMBOL",
                    pvalues = FALSE,
                    nperm = 100,
                    ...) {
    consistencyCheck(nametype, "ECMSign")

    sign_df <- ECM_Chakravarthy
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)

    sign_up <- sign_df[grep("ECMup", sign_df$class), ]
    sign_down <- sign_df[grep("ECMdown", sign_df$class), ]

    percentageOfGenesUsed("ECMSign", datasetm, sign_up$SYMBOL, "up")
    percentageOfGenesUsed("ECMSign", datasetm, sign_down$SYMBOL, "down")

    gene_sets <- list(
        ECM_Chakravarthy_up = sign_up$SYMBOL,
        ECM_Chakravarthy_down = sign_down$SYMBOL
    )

    dots <- list(...)

    args <- matchArguments(
        dots,
        list(
            expr = datasetm,
            gset.idx.list = gene_sets,
            method = "ssgsea",
            kcdf = "Poisson",
            ssgsea.norm = FALSE,
            verbose = FALSE
        )
    )

    gsva_count <- do.call(gsva, args)

    if (pvalues) {
        gsva_pval <- GSVAPvalues(
            expr = datasetm,
            gset.idx.list = gene_sets,
            gsvaResult = gsva_matrix,
            nperm = nperm,
            args = args
        )
        gsva_matrix <- rbind(gsva_matrix, gsva_pval)
    }

    return(returnAsInput(
        userdata = dataset,
        result = gsva_count,
        SignName = "",
        datasetm
    ))
}


#' Homologous Recombination Deficiency Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the HRDS scores is
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @importFrom matrixStats rowMedians
#' @importFrom stats t.test
#'
#' @examples
#' HRDSSign(dataset = subovse)
#'
#' @export
HRDSSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "HRDSSign")

    sign_df <- HRDS_Lu
    sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

    datasetm <- getMatrix(dataset)
    datasetm_n <- datasetm - rowMedians(datasetm)

    percentageOfGenesUsed("HRDSSign", datasetm, sign_df$SYMBOL)

    HRDS_P <- datasetm_n[intersect(
        row.names(datasetm_n),
        sign_df[sign_df$coeff == 1, ]$SYMBOL
    ), ]
    HRDS_N <- datasetm_n[intersect(
        row.names(datasetm_n),
        sign_df[sign_df$coeff == -1, ]$SYMBOL
    ), ]

    score <- unlist(lapply(seq_len(ncol(datasetm_n)), function(x) {
        tmp <- t.test(HRDS_P[, x], HRDS_N[, x], alternative = "two.sided")
        tmp[["statistic"]]
    }))

    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = "HRDS_Lu",
        datasetm
    ))
}


#' Adult Intestinal Stem Cell Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the ISC scores will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' ISCSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
ISCSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "microarray") {
        consistencyCheck(nametype, "ISCSign")

        sign_df <- ISC_MerlosSuarez
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)
        sign_list <- split(sign_df$SYMBOL, sign_df$class)

        datasetm <- getMatrix(dataset)
        datasetm_n <-
            if (inputType == "rnaseq") {
                log2(datasetm + 0.01)
            } else {
                datasetm
            }
        datasetm_n <- datasetm_n - rowMeans(datasetm_n)

        percentageOfGenesUsed(
            "ISCSign", datasetm_n, sign_list$`ISCEphB2`, "ISCEphB2"
        )
        percentageOfGenesUsed(
            "ISCSign", datasetm_n, sign_list$`LateTA`, "LateTA"
        )
        percentageOfGenesUsed(
            "ISCSign", datasetm_n, sign_list$`ISCLgr5`, "ISCLgr5"
        )
        percentageOfGenesUsed(
            "ISCSign", datasetm_n, sign_list$`Prolif`, "Prolif"
        )

        scores <- sapply(sign_list, function(x) {
            colMeans(datasetm_n[intersect(
                x, row.names(datasetm_n - colMeans(datasetm_n))
            ), ])
        })
        colnames(scores) <-
            paste0("ISC_MerlosSuarez_", colnames(scores))

        return(returnAsInput(
            userdata = dataset,
            result = t(scores),
            SignName = "",
            datasetm
        ))
    }

#' VEGF Signature
#'
#' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the VEGF score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' VEGFSign(dataset = subovse)
#'
#' @export
VEGFSign <- function(dataset, nametype = "SYMBOL") {
    consistencyCheck(nametype, "VEGFSign")

    datasetm <- getMatrix(dataset)
    datasetm_n <- log2(datasetm + 0.01)

    score <- statScore(
        VEGF_Hu$SYMBOL, datasetm_n,
        nametype, "mean", "VEGFSign"
    )

    return(returnAsInput(
        userdata = dataset,
        result = score,
        SignName = "VEGF_Hu",
        datasetm
    ))
}


#' DNA Repair Signature
#'
##' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the Angiogenesis score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' DNArepSign(dataset = subovse, inputType = "rnaseq")
#'
#' @export
DNArepSign <-
    function(dataset,
             nametype = "SYMBOL",
             inputType = "microarray") {
        consistencyCheck(nametype, "DNArepSign")

        sign_df <- DNArep_Kang
        sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)

        datasetm <- getMatrix(dataset)

        percentageOfGenesUsed("DNArepSign", datasetm, sign_df$SYMBOL)

        datasetm_n <- if (inputType == "rnaseq") {
            log2(datasetm + 0.01)
        } else {
            datasetm
        }
        datasetm_n <-
            datasetm_n[row.names(datasetm_n) %in% sign_df$SYMBOL, ]
        datasetm_n <- scale(t(datasetm_n), center = TRUE, scale = FALSE)

        medianexp <- apply(datasetm_n, 2, median)

        genes_h <- intersect(
            colnames(datasetm_n), sign_df[sign_df$class == "high", ]$SYMBOL
        )
        genes_l <- intersect(
            colnames(datasetm_n), sign_df[sign_df$class == "low", ]$SYMBOL
        )

        score <- rowSums(datasetm_n[, genes_h] > medianexp[genes_h]) +
            rowSums(datasetm_n[, genes_l] < medianexp[genes_l])

        return(
            returnAsInput(
                userdata = dataset,
                result = score,
                SignName = "DNArep_Kang",
                datasetm
            )
        )
    }


#' IPSOV Signature
#'
##' @inherit EMTSign description
#' @inheritParams EMTSign
#'
#' @return A SummarizedExperiment object in which the IPSOV score will be
#' added in the \code{\link[SummarizedExperiment]{colData}} section.
#'
#' @examples
#' IPSOVSign(dataset = subovse)
#'
#' @export
IPSOVSign <- function(dataset,
                      nametype = "SYMBOL",
                      pvalues = FALSE,
                      nperm = 100,
                      ...) {
    consistencyCheck(nametype, "IPSOVSign")

    # datasetm <- getMatrix(dataset)
    # datasetm_n <- scale(datasetm, center = TRUE, scale = TRUE)
    #
    # sign_df <- IPSOVdata
    # sign_df$SYMBOL <- geneIDtrans(nametype, sign_df$SYMBOL)
    # sign_df <- sign_df[sign_df$SYMBOL %in% rownames(datasetm),]
    # sign_list <- split(sign_df$SYMBOL, sign_df$class)
    #
    # percentageOfGenesUsed("IPSOVSign", datasetm, sign_df$SYMBOL)
    #
    # dots <- list(...)
    # args <- matchArguments(dots, list(
    #     expr = datasetm_n, gset.idx.list = sign_list,
    #     method = "ssgsea", kcdf = "Poisson",
    #     ssgsea.norm = FALSE, verbose = FALSE))
    # gsva_matrix <- do.call(gsva, args)
    #
    # sign_class <- unique(sign_df[,2:3])
    # score <- coeffScore(sign_class, gsva_matrix, "IPSOVSign")
    #
    # return(returnAsInput(
    #     userdata = dataset, result = score, SignName = "IPSOV", datasetm))

    return(dataset)
}
