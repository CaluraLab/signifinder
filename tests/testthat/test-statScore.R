library(signifinder)
library(testthat)

suppressPackageStartupMessages(library(SummarizedExperiment))

test_that("TLSSign works", {
    rmatrix <- fakeData("TLS")
    myres <- TLSSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("TLS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"TLS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"TLS"], "double")
    expect_message(TLSSign(rmatrix), "100")
})

test_that("ExpandedImmuneSign works", {
    rmatrix <- fakeData("expandedImmune_Ayers")
    myres <- expandedImmuneSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("expandedImmune_Ayers" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"expandedImmune_Ayers"], ncol(assay(myres)))
    expect_type(colData(myres)[,"expandedImmune_Ayers"], "double")
    expect_message(expandedImmuneSign(rmatrix), "100")
})

test_that("IFNSign works", {
    rmatrix <- fakeData("IFN_Ayers")
    myres <- IFNSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IFN_Ayers" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IFN_Ayers"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IFN_Ayers"], "double")
    expect_message(IFNSign(rmatrix), "100")
})

test_that("ImmuneCytSign based on Rooney's work", {
    rmatrix <- fakeData("ImmuneCyt_Rooney")
    myres <- ImmuneCytSign(rmatrix, author = "Rooney")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCyt_Rooney" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCyt_Rooney"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCyt_Rooney"], "double")
    expect_message(ImmuneCytSign(rmatrix, author = "Rooney"), "100")
})

test_that("MitoticIndexSign works", {
    rmatrix <- fakeData("MitoticIndex")
    myres <- mitoticIndexSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("MitoticIndex" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"MitoticIndex"], ncol(assay(myres)))
    expect_type(colData(myres)[,"MitoticIndex"], "double")
    expect_message(mitoticIndexSign(rmatrix), "100")
})

test_that("MatrisomeSign works", {
    rmatrix <- fakeData("Matrisome_Yuzhalin")
    myres <- matrisomeSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Matrisome_Yuzhalin" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Matrisome_Yuzhalin"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Matrisome_Yuzhalin"], "integer")
    expect_message(matrisomeSign(rmatrix), "100")
})

test_that("immunoScoreSign based on Roh's work", {
    rmatrix <- fakeData("ImmunoScore_Roh")
    myres <- immunoScoreSign(rmatrix, author = "Roh")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmunoScore_Roh" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmunoScore_Roh"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmunoScore_Roh"], "double")
    expect_message(immunoScoreSign(rmatrix, author = "Roh"), "100")
})

test_that("CINSign works", {
    n <- c("CIN25", "CIN70")
    pname <- sample(n, 1)
    rmatrix <- fakeData(pname)
    myres <- CINSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(all(c("CIN25", "CIN70") %in% colnames(colData(myres))))
    expect_length(colData(myres)[,"CIN25"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CIN70"], "double")
    expect_message(CINSign(rmatrix), "100")
})

test_that("hypoxiaSign works", {
    rmatrix <- fakeData("Hypoxia_Buffa")
    myres <- hypoxiaSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Hypoxia_Buffa" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Hypoxia_Buffa"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Hypoxia_Buffa"], "double")
    expect_message(hypoxiaSign(rmatrix), "100")
})

test_that("CellCycleSign based on Lundberg's work", {
    rmatrix <- fakeData("CellCycle_Lundberg")
    myres <- CellCycleSign(rmatrix, author = "Lundberg")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CellCycle_Lundberg" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CellCycle_Lundberg"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CellCycle_Lundberg"], "integer")
    expect_message(CellCycleSign(rmatrix, author = "Lundberg"), "100")
})

test_that("CellCycleSign based on Davoli's work", {
    rmatrix <- fakeData("CellCycle_Davoli")
    myres <- CellCycleSign(rmatrix, author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CellCycle_Davoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CellCycle_Davoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CellCycle_Davoli"], "double")
    expect_message(CellCycleSign(rmatrix, author = "Davoli"), "100")
})

test_that("VEGFSign works", {
    rmatrix <- fakeData("VEGF_Hu")
    myres <- VEGFSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("VEGF_Hu" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"VEGF_Hu"], ncol(assay(myres)))
    expect_type(colData(myres)[,"VEGF_Hu"], "double")
    expect_message(VEGFSign(rmatrix), "100")
})

test_that("angioSign works", {
    rmatrix <- fakeData("Angiogenesis")
    myres <- angioSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Angiogenesis" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Angiogenesis"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Angiogenesis"], "double")
    expect_message(angioSign(rmatrix), "100")
})

test_that("ImmuneCytSign based on Dabvoli's work", {
    rmatrix <- fakeData("ImmuneCyt_Davoli")
    myres <- ImmuneCytSign(rmatrix, author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCyt_Davoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCyt_Davoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCyt_Davoli"], "double")
    expect_message(ImmuneCytSign(rmatrix, author = "Davoli"), "100")
})

