
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

plot_conditional_hists <- function (midp, seg, ...) {
    layout(1:3)
    par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
    # lims <- quantile(abs(seg), 0.98)
    lims <- max(abs(seg))
    for (q in c(.1, .5, .9)) {
        ut <- (abs(rank(midp)/length(midp) - q) < 0.025)
        outl <- (seg > -lims) & (seg < lims)
        hist(seg[ut & outl], breaks=100, main=sprintf("midparent %0.0f%%", 100*q),
             # xlab="segregation noise",
             xlab='',
             xlim=c(-1,1) * lims,
            ...
        )
        legend("bottomleft", title=sprintf("tail: %.2f%% < %.0f", 100 * mean(seg[ut] < -lims), -lims), legend='', bty='n')
        legend("bottomright", title=sprintf("tail: %.2f%% > %.0f", 100 * mean(seg[ut] > lims), lims), legend='', bty='n')
    }
    mtext("segregation noise", 1, line=2.5, cex=0.75)
}

plot_conditional_ratios <- function (midp, seg, do_legend=TRUE, title=NA, ylim=c(0,3), ...) {
    pbreaks <- pmax(0.01, pmin(0.99, seq(0, 1, length.out=51)))
    qbreaks <- c(-Inf, quantile(seg, pbreaks), Inf)
    total_xh <- hist(seg, breaks=qbreaks, plot=FALSE)
    uk <- seq(2, length(qbreaks)-1)
    par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
    plot(0, 0, type='n', xlim=range(qbreaks, finite=TRUE), ylim=ylim,
         xlab="", ylab="enrichment", ...)
    cols = c(2, 1, 3)
    qvals <- c(.1, .5, .9)
    for (k in seq_along(qvals)) {
        q <- qvals[k]
        ut <- (abs(rank(midp)/length(midp) - q) < 0.025)
        xh <- hist(seg[ut], breaks=qbreaks, plot=FALSE)
        lines(total_xh$mids[uk], xh$density[uk] / total_xh$density[uk], col=cols[k])
    }
    mtext("segregation noise", 1, line=2.5, cex=0.75)
    if (do_legend) legend("top", lty=1, lwd=2, col=cols, legend=sprintf("midparent q=%.0f%%", 100*qvals), title=title)
}
