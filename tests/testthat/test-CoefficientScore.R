library(signifinder)
suppressPackageStartupMessages(library(SummarizedExperiment))

test_that("PyroptosisSign works", {
    pyrnames <- c(
        "Pyroptosis_Ye", "Pyroptosis_Shao", "Pyroptosis_Lin", "Pyroptosis_Li")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    myres <- pyroptosisSign(rmatrix, author = substring(pname, 12))
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(pyroptosisSign(
        rmatrix, author = substring(pname, 12)), "100")
})

test_that("FerroptosysSign work", {
    ferrnames <- c("Ferroptosis_Liang", "Ferroptosis_Li",
                   "Ferroptosis_Liu", "Ferroptosis_Ye")
    fname <- sample(ferrnames, 1)
    rmatrix <- signifinder:::fakeData(fname)
    myres <- ferroptosisSign(rmatrix, author = substring(fname, 13))
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(fname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,fname], ncol(assay(myres)))
    expect_type(colData(myres)[,fname], "double")
    expect_message(ferroptosisSign(
        rmatrix, author = substring(fname, 13)), "100")
})

test_that("LipidMetabolism work", {
    rmatrix <- fakeData("LipidMetabolism")
    myres <- lipidMetabolismSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("LipidMetabolism" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"LipidMetabolism"], ncol(assay(myres)))
    expect_type(colData(myres)[,"LipidMetabolism"], "double")
    expect_message(lipidMetabolismSign(rmatrix), "100")
})

test_that("CD49BSC work", {
    rmatrix <- fakeData("CD49BSC")
    myres <- CD49BSCSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CD49BSC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CD49BSC"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CD49BSC"], "double")
    expect_message(CD49BSCSign(rmatrix), "100")
})

test_that("glycolysisSign works", {
    pyrnames <- c("GlycolysisJiang", "GlycolysisZhangL", "GlycolysisLiu",
                  "GlycolysisYu", "GlycolysisXu", "GlycolysisZhangC")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)

    author <- if(pname == "GlycolysisZhangL"|pname == "GlycolysisZhangC"){
        author <-"Zhang"} else substring(pname, 11)

    myres <- glycolysisSign(rmatrix, author = author)

    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(glycolysisSign(
        rmatrix, author = author), "100")
})

test_that("autophagySign works", {
    pyrnames <- c("AutophagyZhang","AutophagyYue", "AutophagyXu",
                  "AutophagyWang", "AutophagyChenM", "AutophagyHu",
                  "AutophagyHou", "AutophagyFei", "AutophagyFang",
                  "AutophagyChenH")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)

    author <- if(pname == "AutophagyChenM"|pname == "AutophagyChenH"){
        "Chen"} else substring(pname, 10)

    myres <- autophagySign(rmatrix, author = author)

    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    expect_message(autophagySign(
        rmatrix, author = author), "100")
})

