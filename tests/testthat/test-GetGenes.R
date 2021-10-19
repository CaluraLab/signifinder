library(signifinder)

test_that("GetGenes works properly", {
    sname <- sample(SignatureNames, 1)
    expect_equal(ncol(GetGenes(sname)), 2)
    expect_equal(length(unique(GetGenes(sname)[,2])), 1)
    expect_equal(unique(GetGenes(sname)[,2]), sname)
})
