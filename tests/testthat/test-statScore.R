library(signifinder)
library(SummarizedExperiment)
library(labstatR)

test_that("TLSSign works", {
    rmatrix <- fakeData("TLS")
    myres <- TLSSign(rmatrix, nametype = "SYMBOL", tumorTissue = "skin")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("TLS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"TLS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"TLS"], "double")
})

test_that("ExpandedImmuneSign works", {
    rmatrix <- fakeData("ExpandedImmune")
    myres <- expandedImmuneSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ExpandedImmune" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ExpandedImmune"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ExpandedImmune"], "double")
})

test_that("IFNSign works", {
    rmatrix <- fakeData("IFN")
    myres <- IFNSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IFN" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IFN"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IFN"], "double")
})

test_that("CYTSign works", {
    rmatrix <- fakeData("CYT")
    myres <- CYTSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CYT" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CYT"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CYT"], "double")
})

test_that("MitoticIndexSign works", {
    rmatrix <- fakeData("MitoticIndex")
    myres <- mitoticIndexSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("MitoticIndex" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"MitoticIndex"], ncol(assay(myres)))
    expect_type(colData(myres)[,"MitoticIndex"], "double")
})

test_that("MatrisomeSign works", {
    rmatrix <- fakeData("Matrisome")
    myres <- matrisomeSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Matrisome" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Matrisome"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Matrisome"], "integer")
})

test_that("immunoScoreSign works", {
    rmatrix <- fakeData("ImmunoScoreRoh")
    myres <- immunoScoreSign(rmatrix, author = "Roh")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmunoScoreRoh" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmunoScoreRoh"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmunoScoreRoh"], "double")
})


