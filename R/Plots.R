

#' Scatterplot for a single signature
#'
#' Given signatures' scores, it returns a scatterplot of samples' scores and a
#' barplot of the density distributions of samples' scores.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param whichSign character string saying the signature to plot.
#' @param statistics character string saying the statistics to be plotted in the
#' graph. Either one of "mean", "median" or "quantiles".
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median quantile
#'
#' @examples
#' data(ovse)
#' oneSignPlot(data = ovse, whichSign = "Ferroptosis_Ye")
#'
#' @export
oneSignPlot <- function(data, whichSign, statistics = NULL) {

    ..density.. <- NULL

    if (length(whichSign) > 1) {
        stop("you must provide only one signature for this plot")
    }

    .signatureNameCheck(data, whichSign)

    if (!(is.null(statistics))) {
        if (!(statistics %in% c("mean", "median", "quantiles"))) {
            stop("statistics must be one of: mean, median and quantiles.")
        }
    }

    signval <- sort(data.frame(colData(data))[, whichSign])

    g1 <- ggplot() +
        geom_point(mapping = aes(signval, seq_along(signval)), size = 1) +
        labs(x = whichSign, y = "Sample") +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = "grey50")
        )
    g2 <- ggplot() +
        geom_histogram(
            mapping = aes(signval, ..density..),
            colour = "white",
            fill = "skyblue2"
        ) +
        geom_density(mapping = aes(signval), size = 1.5) +
        labs(x = whichSign, y = "Density") +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = "grey50")
        )

    if (!is.null(statistics)) {
        if (statistics == "mean") {
            g1 <- g1 + geom_vline(
                mapping = aes(xintercept = mean(signval)),
                col = "red",
                size = 1.2,
                linetype = 2
            )
            g2 <- g2 + geom_vline(
                mapping = aes(xintercept = mean(signval)),
                col = "red", size = 1.2, linetype = 2)
        } else if (statistics == "median") {
            g1 <- g1 + geom_vline(
                mapping = aes(xintercept = median(signval)),
                col = "red",
                size = 1.2,
                linetype = 2
            )
            g2 <- g2 + geom_vline(
                mapping = aes(xintercept = median(signval)),
                col = "red", size = 1.2, linetype = 2)
        } else if (statistics == "quantiles") {
            g1 <- g1 +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.25)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                ) +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.50)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                ) +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.75)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                )
            g2 <- g2 +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.25)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                ) +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.50)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                ) +
                geom_vline(
                    mapping = aes(xintercept = quantile(signval, 0.75)),
                    col = "red",
                    size = 1.2,
                    linetype = 2
                )
        }
    }
    g1 <- g1 + annotate(
        "text", label = statistics, x = quantile(signval, 0.10),
        y = length(signval), size = 4, colour = "red")

    return(g1 + g2)
}


#' Genes' Signatures' Heatmap
#'
#' Given one or multiple signatures, the function returns a heatmap of the
#' expression values of the genes included in each of them.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param whichSign character vector saying the signatures to plot.
#' @param logCount logical. If TRUE it shows logarithms of expression values.
#' @param splitBySign logical. If TRUE it splits rows by signatures.
#' @param sampleAnnot vector containing samples' annotations.
#' @param splitBySampleAnnot logical. If TRUE it splits columns by samples'
#' annotations.
#' @param ... other parameters specific of the function
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation '%v%' HeatmapAnnotation
#' @importFrom magrittr '%>%'
#' @importFrom dplyr group_by summarise_all
#' @importFrom SummarizedExperiment colData
#' @importFrom grid gpar
#'
#' @examples
#' data(ovse)
#' geneHeatmapSignPlot(data = ovse, whichSign = "Ferroptosis_Ye")
#'
#' @export
geneHeatmapSignPlot <- function(data,
                                whichSign,
                                logCount = FALSE,
                                splitBySign = FALSE,
                                sampleAnnot = NULL,
                                splitBySampleAnnot = FALSE,
                                ...) {
    .signatureNameCheck(data, whichSign)

    dataset <- .getMatrix(data)

    if (!is.null(sampleAnnot)) {
        if (length(sampleAnnot) != ncol(dataset)) {
            stop("sampleAnnot length is different than samples dimension")
        }
    } else {
        if (splitBySampleAnnot) {
            stop("splitBySampleAnnot can be TRUE",
                 " only if sampleAnnot is provided")
        }
    }

    signval <- colData(data)[, whichSign]
    if (length(whichSign) == 1) {
        signval <- matrix(signval,
            nrow = 1,
            dimnames = list(whichSign, colnames(dataset))
        )
    } else {
        signval <- vapply(signval, .range01, double(nrow(signval)))
        signval <- as.matrix(t(signval))
    }

    Gene <- NULL
    geneTable <- as.data.frame(do.call(rbind, lapply(whichSign, .GetGenes))) %>%
        group_by(Gene) %>%
        summarise_all(paste, collapse = ",")
    signatureGenes <- geneTable$Gene

    filtdataset <- as.matrix(dataset[row.names(dataset) %in% signatureGenes, ])

    dots <- list(...)
    htargs <- .matchArguments(
        dots,
        list(
            name = "gene\nexpression",
            show_column_names = FALSE,
            col = mycol,
            row_names_gp = gpar(fontsize = 6)
        )
    )

    if (logCount) {
        htargs$matrix <- log2(filtdataset + 1)
    } else {
        htargs$matrix <- filtdataset
    }

    if (length(whichSign) != 1) {
        signAnnot <- geneTable$Signature[
            geneTable$Gene %in% rownames(filtdataset)]
        if (splitBySign) {
            htargs$row_split <- signAnnot
        } else {
            ha <- rowAnnotation(signature = signAnnot)
            htargs$right_annotation <- ha
        }
    }

    if (splitBySampleAnnot & is.character(sampleAnnot)) {
        htargs$column_split <- sampleAnnot
        ht <- Heatmap(
            signval,
            name = "score",
            col = mycol1,
            column_split = sampleAnnot
        )
    } else {
        if (!is.null(sampleAnnot)) {
            hatop <- HeatmapAnnotation(sampleAnnot = sampleAnnot)
            htargs$top_annotation <- hatop
        }
        ht <- Heatmap(signval, name = "score", col = mycol1)
    }

    ht2 <- do.call(Heatmap, htargs)
    g <- ht %v% ht2

    return(g)
}


#' Global Heatmap of Signatures' scores.
#'
#' Given one or multiple signatures, the function returns a heatmap of scores.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param whichSign character vector saying the signatures to plot. If not
#' specified, all the signatures inside data will be plotted.
#' @param clusterBySign character vector saying one or more signatures to use to
#' cluster columns.
#' @param sampleAnnot vector containing samples' annotations.
#' @param signAnnot character vector of signature's annotations. One or more
#' between: "signature", "topic", "tumor", "tissue".
#' @param splitBySampleAnnot logical. If TRUE it splits columns by samples'
#' annotations.
#' @param ... other parameters specific of the function
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @importFrom ComplexHeatmap Heatmap '%v%'
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' data(ovse)
#' heatmapSignPlot(data = ovse)
#'
#' @export
heatmapSignPlot <- function(
        data, whichSign = NULL, clusterBySign = NULL, sampleAnnot = NULL,
        signAnnot = NULL, splitBySampleAnnot = FALSE, ...) {
    if (!is.null(whichSign)) { .signatureNameCheck(data, whichSign) }
    if (!is.null(clusterBySign)) { .signatureNameCheck(data, clusterBySign) }

    if (sum(colnames(colData(data)) %in% SignatureNames) > 0) {
        data <- colData(data)[, colnames(colData(data)) %in% SignatureNames]
    } else { stop("There are no signatures in data") }

    if (!is.null(sampleAnnot)) {
        if (length(sampleAnnot) != nrow(data)) {
            stop("sampleAnnot length is different than samples dimension")}
    } else {
        if (splitBySampleAnnot) { stop(
            "splitBySampleAnnot can be TRUE only",
            " if sampleAnnot is provided")}}

    if (!is.null(signAnnot)) {
        if (!(signAnnot %in% c("signature", "topic", "tumor", "tissue"))) {
            stop("signAnnot must be one of: signature, topic, tumor, tissue.")
        }
    }

    if (!is.null(whichSign)) {
        data <- data[, intersect(c(whichSign, clusterBySign), colnames(data))]
    }
    keepnames <- rownames(data)

    data <- vapply(data, .range01, double(nrow(data)))
    row.names(data) <- keepnames
    data <- as.matrix(t(data))

    dots <- list(...)
    htargs <- .matchArguments(
        dots,
        list(
            name = "scaled\nscore",
            show_column_names = FALSE,
            col = mycol
        )
    )

    if (!is.null(sampleAnnot)) {
        if (splitBySampleAnnot) {
            htargs$column_split <- sampleAnnot
        } else {
            hatop <- HeatmapAnnotation(sampleAnnot = sampleAnnot)
            htargs$top_annotation <- hatop
        }
    }

    if (is.null(clusterBySign)) {
        if (!is.null(signAnnot)) {
            whichRow <- vapply(
                rownames(data), grep,
                x = signatureTable$scoreLabel,
                FUN.VALUE = integer(1))
            df <- as.data.frame(signatureTable[whichRow, signAnnot])
            colnames(df) <- signAnnot
            ha <- rowAnnotation(df = df)
            htargs$right_annotation <- ha
        }
        htargs$matrix <- data
        g <- do.call(Heatmap, htargs)
    } else {
        n <- which(rownames(data) %in% clusterBySign)
        fm <- as.matrix(data.frame(data)[n, ])
        sm <- as.matrix(data.frame(data)[-n, ])
        if (!is.null(signAnnot)) {
            whichRow <- vapply(
                rownames(sm), grep,
                x = signatureTable$scoreLabel,
                FUN.VALUE = integer(1))
            df <- as.data.frame(signatureTable[whichRow, signAnnot])
            colnames(df) <- signAnnot
            ha <- rowAnnotation(df = df)
            htargs$right_annotation <- ha
        }
        htargs$matrix <- sm
        if (splitBySampleAnnot) {
            ht <- Heatmap(
                fm,
                name = "clustered\nscore",
                col = mycol1,
                column_split = sampleAnnot
            )
        } else {
            ht <- Heatmap(fm, name = "clustered\nscore", col = mycol1)
        }
        g <- ht %v% do.call(Heatmap, htargs)
    }
    return(g)
}


#' Correlation Plot
#'
#' Given multiple signatures, the function plots signatures correlations.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param whichSign character vector saying the signatures to plot. If not
#' specified, all the signatures inside data will be plotted.
#' @param sampleAnnot character vector containing samples' annotations.
#' @param selectByAnnot character string saying the subgroup from `sampleAnnot`
#' used to compute the correlation plot.
#'
#' @return An object of class "openair".
#'
#' @importFrom openair corPlot
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' data(ovse)
#' correlationSignPlot(data = ovse)
#'
#' @export
correlationSignPlot <- function(
        data, whichSign = NULL, sampleAnnot = NULL, selectByAnnot = NULL) {
    if (!is.null(whichSign)) { .signatureNameCheck(data, whichSign) }

    tmp <- colData(data)

    if (sum(colnames(tmp) %in% SignatureNames) > 0) {
        if (is.null(whichSign)) {
            signs <- intersect(SignatureNames, colnames(tmp))
        } else { signs <- Reduce(
            intersect, list(whichSign, SignatureNames, colnames(tmp)))
        }
    } else { stop("There are no signatures in data") }

    tmp <- tmp[, signs]

    if (!is.null(sampleAnnot)) {
        if (length(sampleAnnot) != nrow(tmp)) { stop(
            "sampleAnnot length is different than samples dimension")}
        if (!is.null(selectByAnnot)) {
            if (!(selectByAnnot %in% sampleAnnot)) { stop(
                "selectByAnnot is not present in sampleAnnot")}
        } else { stop(
            "sampleAnnot can be used only if",
            " selectByAnnot is also provided")}
    } else {
        if (!is.null(selectByAnnot)) { stop(
            "selectByAnnot can be used only",
            " if sampleAnnot is also provided")}}

    if (!is.null(sampleAnnot)) {
        if (!is.null(selectByAnnot)) {
            tmp <- tmp[sampleAnnot == selectByAnnot, ] }}

    SignMatrix <- vapply(tmp, .range01, double(nrow(tmp)))

    g <- corPlot(
        as.data.frame(SignMatrix), cluster = TRUE,
        dendrogram = TRUE, lower = TRUE, fontsize = 6)

    return(g)
}


#' Survival Plot
#'
#' Given a signature and samples survival data, the function plots survival
#' curves for that signature.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param survData a dataframe with samples on rows and two columns. The first
#' column holds survival data of time, indicating the follow up times; the
#' second holds data of the survival status, normally 0=alive and 1=dead. For
#' further details check \code{\link[survival]{Surv}} function.
#' @param whichSign character string saying the signature to plot.
#' @param cutpoint a character string (one of: "median", "mean" and "optimal")
#' or a numeric value, which divide samples between high scores and low scores.
#' The function computes the threshold with the method indicated or employs the
#' values directly supplied by the user. Based on that number, it divides
#' samples.
#' @param sampleAnnot a categorical vector containing samples' annotations
#' named with samples names equal to the row names used in `survData`.
#' @param selectByAnnot character string saying the subgroup from `sampleAnnot`
#' used to compute the survival analysis.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median
#' @importFrom maxstat maxstat.test
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#'
#' @examples
#' data(ovse)
#' mysurvData <- cbind(ovse$os, ovse$status)
#' rownames(mysurvData) <- rownames(SummarizedExperiment::colData(ovse))
#' survivalSignPlot(
#'     data = ovse,
#'     survData = mysurvData,
#'     whichSign = "Ferroptosis_Ye"
#' )
#'
#' @export
survivalSignPlot <- function(
        data, survData, whichSign, cutpoint = "mean",
        sampleAnnot = NULL, selectByAnnot = NULL) {
    if (length(whichSign) > 1) {
        stop("you must provide only one signature for this plot")
    }
    .signatureNameCheck(data, whichSign)

    if (ncol(survData) > 2) { stop(
        "survData should contain only two",
        " columns: survival data and status")}

    if (!(is.numeric(cutpoint))) {
        if (!(cutpoint %in% c("mean", "median", "optimal"))) {
            stop(
                "Cutpoint must be either a number or",
                " one of: mean, median and optimal")}}

    tmp <- intersect(rownames(colData(data)), rownames(survData))
    tmp <- as.data.frame(cbind(colData(data)[tmp, ], survData[tmp, ]))

    grp <- rep("high", nrow(tmp))
    names(grp) <- rownames(tmp)
    n <- ncol(tmp)
    colnames(tmp)[c(n - 1, n)] <- c("survival", "status")
    if (cutpoint == "mean") {
        grp[which(tmp[, whichSign] < mean(tmp[, whichSign]))] <- "low"
    } else if (cutpoint == "median") {
        grp[which(tmp[, whichSign] < median(tmp[, whichSign]))] <- "low"
    } else if (cutpoint == "optimal") {
        optval <- maxstat.test(
            Surv(survival, status) ~ tmp[, whichSign],
            data = tmp, smethod = "LogRank", pmethod = "Lau94")
        grp[which(tmp[, whichSign] < optval$estimate)] <- "low"
    } else {
        grp[which(tmp[, whichSign] < cutpoint)] <- "low"
    }

    if ((sum(grp == "low") < length(grp) / 10) |
        (sum(grp == "low") > length(grp) * 9 / 10)) {
        warning(
            "groups size is non homogeneous: ",
            sum(grp == "low"), " low and ",
            sum(grp == "high"), " high")
    }

    if (!is.null(sampleAnnot)) {
        if (length(sampleAnnot) != nrow(tmp)) {
            stop("sampleAnnot length is different than samples dimension")
        }
        if (!is.null(selectByAnnot)) {
            if (!(selectByAnnot %in% sampleAnnot)) {
                stop("selectByAnnot is not present in sampleAnnot")}
        } else { stop(
            "sampleAnnot can be used only if",
            " selectByAnnot is also provided")}
    } else {
        if (!is.null(selectByAnnot)) { stop(
            "selectByAnnot can be used only",
            " if sampleAnnot is also provided")}
    }

    tmp <- cbind(tmp, grp)
    if (!is.null(sampleAnnot)) {
        if (!is.null(selectByAnnot)) {
            tmp <- tmp[sampleAnnot == selectByAnnot, ]}}
    fit <- survfit(Surv(survival, status) ~ grp, data = tmp)
    g <- ggsurvplot(
        fit, data = tmp, risk.table = TRUE, legend.title = whichSign,
        palette = c("red", "blue"), ggtheme = theme_gray(15),
        font.legend = 15, font.tickslab = 15, font.x = 15, font.y = 15,
        risk.table.fontsize = 5, pval = TRUE, surv.median.line = "hv",
        risk.table.col = "strata", tables.height = 0.4)
    return(g)
}


#' Ridgeline Plot
#'
#' Given multiple signatures, the function plots densities scores.
#'
#' @param data an object of type \linkS4class{SummarizedExperiment}. Output of
#' the signatures functions.
#' @param whichSign character vector saying the signatures to plot. If not
#' specified, all the signatures inside data will be plotted.
#' @param groupByAnnot character vector containing samples' annotations.
#' @param selectByAnnot character string saying the subgroup from `groupByAnnot`
#' used to compute the ridgeline plot.
#' @param ... other parameters specific of the functions
#' \code{\link[ggridges]{geom_density_ridges}} and
#' \code{\link[ggridges]{geom_density_ridges_gradient}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges geom_density_ridges_gradient
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' data(ovse)
#' ridgelineSignPlot(data = ovse)
#'
#' @export
ridgelineSignPlot <- function(
        data, whichSign = NULL, groupByAnnot = NULL,
        selectByAnnot = NULL, ...) {
    if (!is.null(whichSign)) {
        .signatureNameCheck(data, whichSign)
    }

    tmp <- colData(data)

    if (sum(colnames(tmp) %in% SignatureNames) > 0) {
        if (is.null(whichSign)) {
            signs <- intersect(SignatureNames, colnames(tmp))
        } else {
            signs <- Reduce(
                intersect,
                list(whichSign, SignatureNames, colnames(tmp))
            )
        }
    } else {
        stop("There are no signatures in data")
    }

    tmp <- as.data.frame(tmp[, signs])
    colnames(tmp) <- signs
    tmp <- data.frame(vapply(tmp, .range01, double(nrow(tmp))))

    if (!is.null(groupByAnnot)) {
        if (length(groupByAnnot) != nrow(tmp)) {
            stop("groupByAnnot length is different than samples dimension")}
        if (!is.null(selectByAnnot)) {
            if (!all(selectByAnnot %in% groupByAnnot)) {
                stop("selectByAnnot is not present in groupByAnnot")}}
    } else {
        if (!is.null(selectByAnnot)) { stop(
            "selectByAnnot can be used only",
            " if groupByAnnot is also provided")}}

    if (is.null(whichSign)) {n <- ncol(tmp)} else {n <- length(whichSign)}

    if (!is.null(selectByAnnot)) {
        tmp <- tmp[groupByAnnot %in% selectByAnnot, ]
        groupByAnnot <- groupByAnnot[groupByAnnot %in% selectByAnnot]
    }

    score <- NULL
    signature <- NULL
    tmp1 <- do.call(rbind, lapply(seq_len(ncol(tmp)), function(x) {
        data.frame(
            score = tmp[, x], signature = colnames(tmp[x]), row.names = NULL)
    }))

    dots <- list(...)
    ridgeargs <- .matchArguments(dots, list(
        alpha = 0.5, bandwidth = 0.05, scale = 1))

    if (is.null(groupByAnnot)) {
        x <- NULL
        g <- ggplot(tmp1, aes(x = score, y = signature, fill = stat(x))) +
            do.call(geom_density_ridges_gradient, ridgeargs) +
            scale_fill_viridis_c(name = "score", option = "A")
    } else {
        ridgeargs$mapping <- aes(fill = rep(groupByAnnot, n))
        g <- ggplot(tmp1, aes(x = score, y = signature)) +
            do.call(geom_density_ridges, ridgeargs) +
            scale_fill_discrete(name = "Group")
    }
        return(g)
}
