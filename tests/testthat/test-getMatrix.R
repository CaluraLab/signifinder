context("Input and output of signature functions")
library(signifinder)

test_that("getMatrix returns matrix", {
    expect_is(getMatrix(matrix(1:20)), "matrix")
    expect_is(getMatrix(data.frame(1:20)), "matrix")
})
