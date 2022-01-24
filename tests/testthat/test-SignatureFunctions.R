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
    myoutput <- capture.output(
        EMTSign(
            rmatrix, nametype = "SYMBOL", tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("EMTSign based on Mak's work", {
    rmatrix  <- fakeData("EMTMak")
    myresMak <- EMTSign(rmatrix, nametype = "SYMBOL",
                    tumorTissue = "pan-tissue", author = "Mak")
    expect_true(is(myresMak, "SummarizedExperiment"))
    expect_true("EMTMak" %in% colnames(colData(myresMak)))
    expect_length(colData(myresMak)[,"EMTMak"], ncol(assay(myresMak)))
    expect_type(colData(myresMak)[,"EMTMak"], "double")
    myoutput <- capture.output(EMTSign(rmatrix, nametype = "SYMBOL",
                                   tumorTissue = "pan-tissue",
                                   author = "Mak"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                      regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("EMTSign based on Cheng's work", {
    rmatrix  <- fakeData("EMTCheng")
    myres <- EMTSign(rmatrix, nametype = "SYMBOL",
                        tumorTissue = "breast", author = "Cheng")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("EMTCheng" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"EMTCheng"], ncol(assay(myres)))
    expect_type(colData(myres)[,"EMTCheng"], "double")
    myoutput <- capture.output(EMTSign(rmatrix, nametype = "SYMBOL",
                                       tumorTissue = "breast",
                                       author = "Cheng"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
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
    myoutput <- capture.output(
        platinumResistanceSign(
            rmatrix, nametype = "SYMBOL", tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("ASCSign work", {
    rmatrix  <- signifinder:::fakeData("ASC")
    myres <- ASCSign(rmatrix, nametype = "SYMBOL", tumorTissue = "epithelial-derived neuroendocrine cancer")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ASC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ASC"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ASC"], "double")
    myoutput <- capture.output(ASCSign(rmatrix, nametype = "SYMBOL",
                                       tumorTissue = "epithelial-derived neuroendocrine cancer"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("chemokineSign work", {
    rmatrix  <- signifinder:::fakeData("Chemokine")
    myres <- chemokineSign(rmatrix, nametype = "SYMBOL", tumorTissue = "pan-tissue")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Chemokine" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Chemokine"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Chemokine"], "double")
    myoutput <- capture.output(chemokineSign(rmatrix, nametype = "SYMBOL",
                                             tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

# test_that("PassONSign work", {
#     rmatrix  <- signifinder:::fakeData("PASS.ON")
#     myres <- PassONSign(rmatrix, nametype = "SYMBOL", tumorTissue = "skin")
#     expect_true(is(myres, "SummarizedExperiment"))
#     expect_true("Pass.ON" %in% colnames(colData(myres)))
#     expect_length(colData(myres)[,"Pass.ON"], ncol(assay(myres)))
#     expect_type(colData(myres)[,"Pass.ON"], "double")
#     myoutput <- capture.output(PassONSign(rmatrix, nametype = "SYMBOL",
#                                              tumorTissue = "skin"))[1]
#     myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
#                           regexpr("%", myoutput)-1)
#     expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
# })


test_that("CISSign work", {
    pyrnames <- c("CISup", "CISdown")
    pname <- sample(pyrnames, 1)
    rmatrix  <- signifinder:::fakeData(pname)
    myres <- CISSign(rmatrix, nametype = "SYMBOL", tumorTissue = "bladder")
    expect_true(is(myres, "SummarizedExperiment"))
    pname <- substring(pname, 1,3)
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    myoutput <- capture.output(CISSign(rmatrix, nametype = "SYMBOL",
                                       tumorTissue = "bladder"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("HRDSSign work", {
    rmatrix  <- signifinder:::fakeData("HRDS")
    myres <- HRDSSign(rmatrix, nametype = "SYMBOL", tumorTissue = "ovary")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("HRDS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"HRDS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"HRDS"], "double")
    myoutput <- capture.output(HRDSSign(rmatrix, nametype = "SYMBOL",
                                             tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("DNArepSign work", {
    rmatrix  <- signifinder:::fakeData("DNArepair")
    myres <- DNArepSign(rmatrix, nametype = "SYMBOL", tumorTissue = "ovary")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("DNArepair" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"DNArepair"], ncol(assay(myres)))
    expect_type(colData(myres)[,"DNArepair"], "integer")
    myoutput <- capture.output(DNArepSign(rmatrix, nametype = "SYMBOL",
                                             tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("IPRES work", {
    # rmatrix  <- fakeData("IPRES")
    myres <- IPRESSign(datasetm, nametype = "SYMBOL", tumorTissue = "skin")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IPRES" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IPRES"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IPRES"], "double")
    myoutput <- capture.output(IPRESSign(datasetm, nametype = "SYMBOL",
                                       tumorTissue = "skin"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})
