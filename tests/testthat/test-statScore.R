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

# test_that("CYTSign works", {
#     rmatrix <- fakeData("CYT")
#     myres <- CYTSign(rmatrix)
#     expect_true(is(myres, "SummarizedExperiment"))
#     expect_true("CYT" %in% colnames(colData(myres)))
#     expect_length(colData(myres)[,"CYT"], ncol(assay(myres)))
#     expect_type(colData(myres)[,"CYT"], "double")
#     myoutput <- capture.output(CYTSign(rmatrix, nametype = "SYMBOL",
#                                         tumorTissue = "pan-tissue"))[1]
#     myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
#                             regexpr("%", myoutput)-1)
#     expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
# })

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

test_that("immunoScoreSign works", {
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
    expect_true("CIN" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CIN"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CIN"], "integer")
    myoutput <- capture.output(CINSign(rmatrix, nametype = "SYMBOL",
            tumorTissue = "breast, lung, brain, lymphatic system"))[1]
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

test_that("CCSSign works", {
    rmatrix <- signifinder:::fakeData("CCSLundberg")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Lundberg")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSLundberg" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSLundberg"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSLundberg"], "integer")
    myoutput <- capture.output(CCSSign(rmatrix, nametype = "SYMBOL"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                            regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})


test_that("CCSSign works", {
    rmatrix <- signifinder:::fakeData("CCSDavoli")
    myres <- CCSSign(rmatrix, tumorTissue = "pan-tissue", author = "Davoli")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CCSDavoli" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CCSDavoli"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CCSDavoli"], "double")
    myoutput <- capture.output(CCSSign(rmatrix, nametype = "SYMBOL"))[1]
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
    myoutput <- capture.output(VEGFSign(rmatrix, nametype = "SYMBOL"))[1]
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
    myoutput <- capture.output(angioSign(rmatrix, nametype = "SYMBOL"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("CyISign works", {
    rmatrix <- signifinder:::fakeData("CytoImmuno")
    myres <- CyISign(rmatrix, tumorTissue = "pan-tissue")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CytoImmuno" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CytoImmuno"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CytoImmuno"], "double")
    myoutput <- capture.output(CyISign(rmatrix, nametype = "SYMBOL"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})
