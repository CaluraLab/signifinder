library(signifinder)
suppressPackageStartupMessages(library(SummarizedExperiment))

test_that("PyroptosisSign works", {
    pyrnames <- c(
        "PyroptosisYe", "PyroptosisShao", "PyroptosisLin", "PyroptosisLi")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    tissue <- signatureTable$tumorTissue[
        signatureTable$functionName=="pyroptosisSign" &
            signatureTable$author==substring(pname, 11)]
    myres <- pyroptosisSign(
        rmatrix, tumorTissue = tissue, author = substring(pname, 11))
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    myoutput <- capture.output(pyroptosisSign(
        rmatrix, nametype = "SYMBOL", tumorTissue = tissue, author = substring(
            pname, 11)))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("FerroptosysSign work", {
    ferrnames <- c("FerroptosisChang", "FerroptosisLiang", "FerroptosisLi",
                   "FerroptosisLiu", "FerroptosisYe", "FerroptosisZhu")
    fname <- sample(ferrnames, 1)
    rmatrix <- fakeData(fname)
    tissue <- signatureTable$tumorTissue[
        signatureTable$functionName=="ferroptosisSign" &
            signatureTable$author==substring(fname, 12)]
    myres <- ferroptosisSign(
        rmatrix, tumorTissue = tissue, author = substring(fname, 12))
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(fname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,fname], ncol(assay(myres)))
    expect_type(colData(myres)[,fname], "double")
    myoutput <- capture.output(ferroptosisSign(
        rmatrix, nametype = "SYMBOL", tumorTissue = tissue,
        author = substring(fname, 12)))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("LipidMetabolism work", {
    rmatrix <- fakeData("LipidMetabolism")
    myres <- lipidMetabolismSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("LipidMetabolism" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"LipidMetabolism"], ncol(assay(myres)))
    expect_type(colData(myres)[,"LipidMetabolism"], "double")
    myoutput <- capture.output(lipidMetabolismSign(rmatrix, nametype = "SYMBOL",
                                                   tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("CD49BSC work", {
    rmatrix <- fakeData("CD49BSC")
    myres <- CD49BSCSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CD49BSC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CD49BSC"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CD49BSC"], "double")
    myoutput <- capture.output(CD49BSCSign(rmatrix, nametype = "SYMBOL",
                                           tumorTissue = "prostate"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("glycolysisSign works", {
    pyrnames <- c("GlycolysisJiang", "GlycolysisZhangL", "GlycolysisLiu",
                  "GlycolysisYu", "GlycolysisXu", "GlycolysisZhangC")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    tissue <- if(pname == "GlycolysisZhangL"){
        tissue <-"lung"} else if (pname == "GlycolysisZhangC"){
            tissue <-"bladder"} else signatureTable$tumorTissue[
        signatureTable$functionName=="glycolysisSign" &
            signatureTable$author==substring(pname, 11)]

    myres <- glycolysisSign(
        rmatrix, tumorTissue = tissue,
        author = if(pname == "GlycolysisZhangL"|pname == "GlycolysisZhangC"){
            author <-"Zhang"} else substring(pname, 11))

    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    myoutput <- capture.output(glycolysisSign(
        rmatrix, nametype = "SYMBOL", tumorTissue = tissue, author = if(
            pname == "GlycolysisZhangL"|pname == "GlycolysisZhangC"){
            author <-"Zhang"} else substring(pname, 11)))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("autophagySign works", {
    pyrnames <- c("AutophagyZhang","AutophagyYue", "AutophagyXu",
                  "AutophagyWang", "AutophagyChenM", "AutophagyHu",
                  "AutophagyHou", "AutophagyFei", "AutophagyFang",
                  "AutophagyChenH")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    tissue <- if(pname == "AutophagyChenM"){
        tissue <-"kidney"} else if (pname == "AutophagyChenH"){
            tissue <-"cervix"} else signatureTable$tumorTissue[
                signatureTable$functionName=="autophagySign" &
                    signatureTable$author==substring(pname, 10)]

    myres <- autophagySign(
        rmatrix, tumorTissue = tissue,
        author = if(pname == "AutophagyChenM"|pname == "AutophagyChenH"){
            author <-"Chen"} else substring(pname, 10))

    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    myoutput <- capture.output(autophagySign(
        rmatrix, nametype = "SYMBOL", tumorTissue = tissue, author = if(
            pname == "AutophagyChenM"|pname == "AutophagyChenH"){
            author <-"Chen"} else substring(pname, 10)))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2,
                          regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

