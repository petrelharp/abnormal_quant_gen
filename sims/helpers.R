
revCumRowSums <- function (x) {
    for (j in ncol(x):2) x[,j] <- rowSums(x[,1:j])
    return(x)
}


plot_tail <- function (x, minp=0.05, ...) {
    x <- abs(x)
    pp <- exp(seq(log(minp), log(.0001), length.out=201))
    qq <- quantile(x, 1 - pp)
    xlm <- lm(log(pp) ~ log(qq))
    plot(qq, pp, log='xy', type='l', ylab="tail probability", ...)
    # abline(coef(xlm), col='red')
    lines(qq, exp(coef(xlm)[1] + coef(xlm)[2] * log(qq)), col='red')
    legend("topright", lty=1, col='red', bty='n',
           legend=sprintf("y = %.2f %.2f x", coef(xlm)[1], coef(xlm)[2])
    )
}
