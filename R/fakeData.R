fakeData <- function(sname){
    g <- GetGenes(sname)
    n <- length(g)*5
    rmatrix <- matrix(rpois(n, 100), ncol = 5)
    colnames(rmatrix) <- g
}
