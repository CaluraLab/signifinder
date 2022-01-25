library(signifinder)
library(testthat)
suppressPackageStartupMessages(library(SummarizedExperiment))

test_that("EMTSign based on Miow's work", {
    pyrnames <- c("Epithelial", "Mesenchymal")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    myres <- EMTSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(EMTSign(rmatrix), "100")
})

test_that("EMTSign based on Mak's work", {
    rmatrix  <- fakeData("EMTMak")
    myresMak <- EMTSign(rmatrix, tumorTissue = "pan-tissue", author = "Mak")
    expect_true(is(myresMak, "SummarizedExperiment"))
    expect_true("EMTMak" %in% colnames(colData(myresMak)))
    expect_length(colData(myresMak)[,"EMTMak"], ncol(assay(myresMak)))
    expect_type(colData(myresMak)[,"EMTMak"], "double")
    expect_message(EMTSign(
        rmatrix, tumorTissue = "pan-tissue", author = "Mak"), "100")
})

test_that("EMTSign based on Cheng's work", {
    rmatrix  <- fakeData("EMTCheng")
    myres <- EMTSign(rmatrix, tumorTissue = "breast", author = "Cheng")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("EMTCheng" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"EMTCheng"], ncol(assay(myres)))
    expect_type(colData(myres)[,"EMTCheng"], "double")
    expect_message(EMTSign(
        rmatrix, tumorTissue = "breast", author = "Cheng"), "100")
})

test_that("platinumResistanceSign work", {
    pyrnames <- c("PlatinumResistanceUp", "PlatinumResistanceDown")
    pname <- sample(pyrnames, 1)
    rmatrix <- signifinder:::fakeData(pname)
    myres <- platinumResistanceSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(platinumResistanceSign(rmatrix), "100")
})

test_that("ASCSign work", {
    rmatrix  <- signifinder:::fakeData("ASC")
    myres <- ASCSign(rmatrix, tumorTissue = "epithelial-derived neuroendocrine cancer")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ASC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ASC"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ASC"], "double")
    expect_message(ASCSign(rmatrix, tumorTissue = "epithelial-derived neuroendocrine cancer"), "100")
})

test_that("chemokineSign work", {
    rmatrix  <- signifinder:::fakeData("Chemokine")
    myres <- chemokineSign(rmatrix, tumorTissue = "pan-tissue")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Chemokine" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Chemokine"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Chemokine"], "double")
    expect_message(chemokineSign(rmatrix, tumorTissue = "pan-tissue"), "100")
})

# test_that("PassONSign work", {
#     rmatrix  <- signifinder:::fakeData("PASS.ON")
#     myres <- PassONSign(rmatrix, tumorTissue = "skin")
#     expect_true(is(myres, "SummarizedExperiment"))
#     expect_true("Pass.ON" %in% colnames(colData(myres)))
#     expect_length(colData(myres)[,"Pass.ON"], ncol(assay(myres)))
#     expect_type(colData(myres)[,"Pass.ON"], "double")
#     expect_message(PassONSign(rmatrix, tumorTissue = "skin"), "100")
# })


test_that("CISSign work", {
    pyrnames <- c("CISup", "CISdown")
    pname <- sample(pyrnames, 1)
    rmatrix  <- signifinder:::fakeData(pname)
    myres <- CISSign(rmatrix, tumorTissue = "bladder")
    expect_true(is(myres, "SummarizedExperiment"))
    pname <- substring(pname, 1,3)
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(CISSign(rmatrix, tumorTissue = "bladder"), "100")
})

test_that("HRDSSign work", {
    rmatrix  <- signifinder:::fakeData("HRDS")
    myres <- HRDSSign(rmatrix, tumorTissue = "ovary")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("HRDS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"HRDS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"HRDS"], "double")
    expect_message(HRDSSign(rmatrix, tumorTissue = "ovary"), "100")
})

test_that("DNArepSign work", {
    rmatrix  <- signifinder:::fakeData("DNArepair")
    myres <- DNArepSign(rmatrix, tumorTissue = "ovary")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("DNArepair" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"DNArepair"], ncol(assay(myres)))
    expect_type(colData(myres)[,"DNArepair"], "integer")
    expect_message(DNArepSign(rmatrix, tumorTissue = "ovary"), "100")
})

# test_that("IPRES work", {
#     rmatrix  <- fakeData("IPRES")
#     myres <- IPRESSign(rmatrix, tumorTissue = "skin")
#     expect_true(is(myres, "SummarizedExperiment"))
#     expect_true("IPRES" %in% colnames(colData(myres)))
#     expect_length(colData(myres)[,"IPRES"], ncol(assay(myres)))
#     expect_type(colData(myres)[,"IPRES"], "double")
#     expect_message(IPRESSign(datasetm, tumorTissue = "skin"), "100")
# })
