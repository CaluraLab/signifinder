library(signifinder)
library(testthat)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(labstatR))

test_that("TLSSign works", {
    rmatrix <- fakeData("TLS")
    myres <- TLSSign(rmatrix, nametype = "SYMBOL", tumorTissue = "skin")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("TLS" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"TLS"], ncol(assay(myres)))
    expect_type(colData(myres)[,"TLS"], "double")
    myoutput <- capture.output(TLSSign(rmatrix, nametype = "SYMBOL",
                                        tumorTissue = "skin"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("ExpandedImmuneSign works", {
    rmatrix <- fakeData("ExpandedImmune")
    myres <- expandedImmuneSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ExpandedImmune" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ExpandedImmune"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ExpandedImmune"], "double")
    myoutput <- capture.output(expandedImmuneSign(rmatrix, nametype = "SYMBOL",
                                            tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("IFNSign works", {
    rmatrix <- fakeData("IFN")
    myres <- IFNSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IFN" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IFN"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IFN"], "double")
    myoutput <- capture.output(IFNSign(rmatrix, nametype = "SYMBOL",
                                        tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("ImmuneCytSign based on Rooney's work", {
    rmatrix <- fakeData("ImmuneCytRooney")
    myres <- ImmuneCytSign(rmatrix,nametype= "SYMBOL",tumorTissue= "pan-tissue",
                           author = "Rooney")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCytRooney" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCytRooney"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCytRooney"], "double")
    myoutput <- capture.output(ImmuneCytSign(rmatrix, nametype = "SYMBOL",
                                             tumorTissue = "pan-tissue",
                                             author = "Rooney"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("MitoticIndexSign works", {
    rmatrix <- fakeData("MitoticIndex")
    myres <- mitoticIndexSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("MitoticIndex" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"MitoticIndex"], ncol(assay(myres)))
    expect_type(colData(myres)[,"MitoticIndex"], "double")
    myoutput <- capture.output(mitoticIndexSign(rmatrix, nametype = "SYMBOL",
                                        tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("MatrisomeSign works", {
    rmatrix <- fakeData("Matrisome")
    myres <- matrisomeSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Matrisome" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Matrisome"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Matrisome"], "integer")
    myoutput <- capture.output(matrisomeSign(rmatrix, nametype = "SYMBOL",
                                                tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("immunoScoreSign based on Roh's work", {
    rmatrix <- fakeData("ImmunoScoreRoh")
    myres <- immunoScoreSign(rmatrix, tumorTissue= "pan-tissue", author= "Roh")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmunoScoreRoh" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmunoScoreRoh"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmunoScoreRoh"], "double")
    myoutput <- capture.output(immunoScoreSign(
        rmatrix, nametype= "SYMBOL", tumorTissue= "pan-tissue",
        author = "Roh"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("CINSign works", {
    rmatrix <- fakeData("CIN")
    myres <- CINSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(all(c("CIN25", "CIN70") %in% colnames(colData(myres))))
    expect_length(colData(myres)[,"CIN25"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CIN70"], "double")
    myoutput <- capture.output(CINSign(rmatrix, nametype = "SYMBOL",
            tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})


test_that("hypoxiaSign works", {
    rmatrix <- fakeData("Hypoxia")
    myres <- hypoxiaSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Hypoxia" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Hypoxia"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Hypoxia"], "double")
    myoutput <- capture.output(hypoxiaSign(rmatrix, nametype = "SYMBOL",
                                               tumorTissue = "pan-tissue"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("CCSSign based on Lundberg's work", {
    rmatrix <- signifinder:::fakeData("CCSLundberg")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Lundberg")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSLundberg" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSLundberg"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSLundberg"], "integer")
    myoutput <- capture.output(CCSSign(rmatrix, tumorTissue = "pan-tissue",
                                       author = "Lundberg"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})


test_that("CCSSign based on Davoli's work", {
    rmatrix <- signifinder:::fakeData("CCSDavoli")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSDavoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSDavoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSDavoli"], "double")
    myoutput <- capture.output(CCSSign(rmatrix, tumorTissue = "pan-tissue",
                                       author = "Davoli"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("VEGFSign works", {
    rmatrix <- signifinder:::fakeData("VEGF")
    myres <- VEGFSign(rmatrix, tumorTissue = "ovary")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("VEGF" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"VEGF"], ncol(assay(myres)))
    expect_type(colData(myres)[,"VEGF"], "double")
    myoutput <- capture.output(VEGFSign(rmatrix, tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("angioSign works", {
    rmatrix <- signifinder:::fakeData("Angiogenesis")
    myres <- angioSign(rmatrix, tumorTissue = "pan-tissue")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Angiogenesis" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Angiogenesis"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Angiogenesis"], "double")
    myoutput <-capture.output(angioSign(rmatrix, tumorTissue = "pan-tissue"))[1]
    myoutput <-substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("ImmuneCytSign based on Dabvoli's work", {
    rmatrix <- fakeData("ImmuneCytDavoli")
    myres <- ImmuneCytSign(rmatrix, nametype="SYMBOL",tumorTissue= "pan-tissue",
                           author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ImmuneCytDavoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"ImmuneCytDavoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"ImmuneCytDavoli"], "double")
    myoutput <- capture.output(ImmuneCytSign(rmatrix, nametype = "SYMBOL",
                                             tumorTissue = "pan-tissue",
                                             author = "Davoli"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

