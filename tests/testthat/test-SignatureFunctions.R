library(signifinder)
suppressPackageStartupMessages(library(SummarizedExperiment))

test_that("EMTSign work", {
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

test_that("EMTSign work", {
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

test_that("EMTSign work", {
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
    rmatrix <- fakeData(pname)
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
