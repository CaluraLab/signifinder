fakeData <- function(sname, input = "rnaseq"){
    g <- GetGenes(sname)
    g <- g[,1]
    n <- length(g)*5
    if(input == "rnaseq"){
        rmatrix <- matrix(rpois(n, 100), ncol = 5)
    } else if (input == "microarray"){
        rmatrix <- matrix(rnorm(n, 100, 50), ncol = 5)
    }
    rownames(rmatrix) <- g
    return(rmatrix)
}
