
revCumRowSums <- function (x) {
    for (j in ncol(x):2) x[,j] <- rowSums(x[,1:j])
    return(x)
}

