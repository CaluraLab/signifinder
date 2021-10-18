context("Output of signature functions")
library(signifinder)

test_that("returnAsInput returns a S4 object", {
    rmatrix <- matrix(rpois(30, 100), ncol = 6)

    expect_true(is.matrix(getMatrix(rmatrix)))
    expect_true(is.matrix(getMatrix(data.frame(rmatrix))))
    expect_equal(typeof(rmatrix), typeof(getMatrix(rmatrix)))
})
