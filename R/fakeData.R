fakeData <- function(sname){
    g <- GetGenes(sname)
    g <- g[,1]
    n <- length(g)*5
    rmatrix <- matrix(rpois(n, 100), ncol = 5)
    rownames(rmatrix) <- g
    return(rmatrix)
}
