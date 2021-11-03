library(signifinder)
library(SummarizedExperiment)

test_that("PyroptosisSign works", {
    pyrnames <- c("PyroptosisYe", "PyroptosisShao", "PyroptosisLin", "PyroptosisLi")
    pname <- sample(pyrnames, 1)
    rmatrix <- fakeData(pname)
    myres <- pyroptosisSign(rmatrix, author = substring(pname, 11))
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[,pname], ncol(assay(myres)))
    expect_type(colData(myres)[,pname], "double")
    myoutput <- capture.output(pyroptosisSign(rmatrix, nametype = "SYMBOL", tumorTissue = "ovary", author = substring(pname, 11)))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2, regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("FerroptosysSign work", {
    rmatrix <- fakeData("Ferroptosis")
    myres <- ferroptosisSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Ferroptosis" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"Ferroptosis"], ncol(assay(myres)))
    expect_type(colData(myres)[,"Ferroptosis"], "double")
    myoutput <- capture.output(ferroptosisSign(rmatrix, nametype = "SYMBOL", tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2, regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("LipidMetabolism work", {
    rmatrix <- fakeData("LipidMetabolism")
    myres <- lipidMetabolismSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("LipidMetabolism" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"LipidMetabolism"], ncol(assay(myres)))
    expect_type(colData(myres)[,"LipidMetabolism"], "double")
    myoutput <- capture.output(lipidMetabolismSign(rmatrix, nametype = "SYMBOL", tumorTissue = "ovary"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2, regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})

test_that("CD49BSC work", {
    rmatrix <- fakeData("CD49BSC")
    myres <- CD49BSCSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CD49BSC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"CD49BSC"], ncol(assay(myres)))
    expect_type(colData(myres)[,"CD49BSC"], "double")
    myoutput <- capture.output(CD49BSCSign(rmatrix, nametype = "SYMBOL", tumorTissue = "prostate"))[1]
    myoutput <- substring(myoutput, regexpr("g ", myoutput)+2, regexpr("%", myoutput)-1)
    expect_true(as.numeric(myoutput)>=0 & as.numeric(myoutput)<=100)
})
