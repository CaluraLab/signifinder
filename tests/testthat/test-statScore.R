library(signifinder)
library(testthat)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(labstatR))

test_that("TLSSign works", {
    rmatrix <- fakeData("TLS")
    myres <- TLSSign(rmatrix, tumorTissue = "skin")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("TLS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"TLS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"TLS"], "double")
    expect_message(TLSSign(rmatrix, tumorTissue = "skin"), "100")
})

test_that("ExpandedImmuneSign works", {
    rmatrix <- fakeData("ExpandedImmune")
    myres <- expandedImmuneSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ExpandedImmune" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ExpandedImmune"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ExpandedImmune"], "double")
    expect_message(expandedImmuneSign(rmatrix), "100")
})

test_that("IFNSign works", {
    rmatrix <- fakeData("IFN")
    myres <- IFNSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IFN" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IFN"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IFN"], "double")
    expect_message(IFNSign(rmatrix), "100")
})

test_that("ImmuneCytSign based on Rooney's work", {
    rmatrix <- fakeData("ImmuneCytRooney")
    myres <- ImmuneCytSign(rmatrix, tumorTissue = "pan-tissue",
                           author = "Rooney")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCytRooney" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCytRooney"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCytRooney"], "double")
    expect_message(ImmuneCytSign(
        rmatrix, tumorTissue = "pan-tissue", author = "Rooney"), "100")
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
    rmatrix <- fakeData("Matrisome")
    myres <- matrisomeSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Matrisome" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Matrisome"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Matrisome"], "integer")
    expect_message(matrisomeSign(rmatrix), "100")
})

test_that("immunoScoreSign based on Roh's work", {
    rmatrix <- fakeData("ImmunoScoreRoh")
    myres <- immunoScoreSign(rmatrix, tumorTissue= "pan-tissue", author= "Roh")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmunoScoreRoh" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmunoScoreRoh"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmunoScoreRoh"], "double")
    expect_message(immunoScoreSign(
        rmatrix, tumorTissue = "pan-tissue", author = "Roh"), "100")
})

test_that("CINSign works", {
    rmatrix <- fakeData("CIN")
    myres <- CINSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(all(c("CIN25", "CIN70") %in% colnames(colData(myres))))
    expect_length(colData(myres)[,"CIN25"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CIN70"], "double")
    expect_message(CINSign(rmatrix), "100")
})

test_that("hypoxiaSign works", {
    rmatrix <- fakeData("Hypoxia")
    myres <- hypoxiaSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Hypoxia" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Hypoxia"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Hypoxia"], "double")
    expect_message(hypoxiaSign(rmatrix), "100")
})

test_that("CCSSign based on Lundberg's work", {
    rmatrix <- signifinder:::fakeData("CCSLundberg")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Lundberg")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSLundberg" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSLundberg"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSLundberg"], "integer")
    expect_message(CCSSign(
        rmatrix, tumorTissue = "pan-tissue", author = "Lundberg"), "100")
})

test_that("CCSSign based on Davoli's work", {
    rmatrix <- signifinder:::fakeData("CCSDavoli")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSDavoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSDavoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSDavoli"], "double")
    expect_message(CCSSign(
        rmatrix, tumorTissue = "pan-tissue", author = "Davoli"), "100")
})

test_that("VEGFSign works", {
    rmatrix <- signifinder:::fakeData("VEGF")
    myres <- VEGFSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("VEGF" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"VEGF"], ncol(assay(myres)))
    expect_type(colData(myres)[,"VEGF"], "double")
    expect_message(VEGFSign(rmatrix), "100")
})

test_that("angioSign works", {
    rmatrix <- signifinder:::fakeData("Angiogenesis")
    myres <- angioSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Angiogenesis" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Angiogenesis"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Angiogenesis"], "double")
    expect_message(angioSign(rmatrix), "100")
})

test_that("ImmuneCytSign based on Dabvoli's work", {
    rmatrix <- fakeData("ImmuneCytDavoli")
    myres <- ImmuneCytSign(rmatrix, tumorTissue= "pan-tissue",
                           author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCytDavoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCytDavoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCytDavoli"], "double")
    expect_message(ImmuneCytSign(
        rmatrix, tumorTissue= "pan-tissue", author = "Davoli"), "100")
})

