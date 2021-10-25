library(signifinder)
library(dplyr)

test_that("signatureTable has not double rows", {
    filtTable <- signatureTable[, c("functionName", "tumorTissue", "author")]
    res <- filtTable %>% filter(duplicated(.))
    expect_equal(nrow(res), 0)
})
