library(signifinder)

test_that("GetGenes works properly", {
    snames <- SignatureNames[!(
        SignatureNames %in% c("MES_consensus", "IMR_consensus",
                              "DIF_consensus", "PRO_consensus"))]
    sname <- sample(snames, 1)
    res <- GetGenes(sname)
    expect_equal(ncol(res), 2)
    expect_equal(length(unique(res[,2])), 1)
    expect_equal(unique(res[,2]), sname)
})
