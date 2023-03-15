# number of individuals
N <- 400

# number of generations
ngens <- 200

# trait dimensionality
trait_dim <- 2

# segregation SD
seg_sd <- 0.5

# mean proportion replaced per step, also equal to dt
alpha <- 0.1

# mean replacement rate depends on 2nd coord
death <- function (x, eta=1) { 1 + (x[,2]/eta)^2 }

xlim <- c(-1,1) 
ylim <- c(-1,1)
plotfun <- function (X, ...) {
    plot(X, pch=20, cex=0.5, col=adjustcolor("black", 0.5), 
         xlab='', ylab='', ...)
    points(t(colMeans(X)), col=adjustcolor('red',0.5), pch=20, cex=2)
    abline(h=0)
}

# initial population mean
init_mean <- c(0,2)

# initial SD
init_sd <- 0.2

# population
X <- sweep(matrix(rnorm(N*trait_dim, sd=init_sd), ncol=trait_dim), 2, init_mean, "+")

t <- 0
plotfun(X, main=sprintf("gen = %0.1f", t*alpha))

for (t in 1:ceiling(ngens/alpha)) {
    dvec <- death(X)
    total_rate <- sum(dvec)
    noff <- min(N, rpois(1, alpha * total_rate))
    i <- sample(nrow(X), size=noff)
    j <- sample(nrow(X), size=noff)
    k <- sample(nrow(X), size=noff, prob=dvec)
    z <- (X[i,,drop=FALSE] + X[j,,drop=FALSE])/2 + matrix(rnorm(trait_dim*noff, sd=seg_sd), ncol=trait_dim)
    X[k,] <- z
    xlim <- range(xlim, X[,1])
    ylim <- range(ylim, X[,2])
    plotfun(X, main=sprintf("gen = %0.1f", t*alpha),
            xlim=xlim, ylim=ylim)
}
