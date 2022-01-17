library(signifinder)
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

test_that("CYTSign works", {
    rmatrix <- fakeData("CYT")
    myres <- CYTSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CYT" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CYT"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CYT"], "double")
    myoutput <- capture.output(CYTSign(rmatrix, nametype = "SYMBOL",
                                        tumorTissue = "pan-tissue"))[1]
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
#
# test_that("CCSSign works", {
#     rmatrix <- fakeData("CCS")
#     myres <- CCSSign(rmatrix)
#     expect_true(is(myres, "SummarizedExperiment"))
#     expect_true("CCS" %in% colnames(colData(myres)))
#     expect_length(colData(myres)[,"CCS"], ncol(assay(myres)))
#     expect_type(colData(myres)[,"CCS"], "integer")
#     myoutput <- capture.output(immunoScoreSign(rmatrix, nametype = "SYMBOL"))[1]
#     myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
#                             regexpr("%", myoutput)-1)
#     expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
# })

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


