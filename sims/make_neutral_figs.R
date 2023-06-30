#!/usr/bin/Rscript

usage <- "
   make_neutral_figs.R (basename)
where basename.repro.tsv and basename.pop.tsv and basename.fix.tsv should be files,
and the 'repro' file should have a first line with parameters,
and remaining lines of the form
  (time) (trait values) (maternal trait values) (paternal trait values)
with one line per individual (NOTE THIS DIFFERS FROM make_figs.R!!);
the 'pop' file should have 'trait' and 'age' columns,
and the 'fix' file should have 'time' and 'num_fixations' columns.
"

library(jsonlite)
library(pracma)

if (!exists("base")) {
    args <- commandArgs(TRUE)
    if (length(args) != 1) {
        stop(usage)
    }

    base <- args[1]
}

infile <- file(paste0(base, ".repro.tsv"), "r")
params <- fromJSON(readLines(infile, 1))
repro <- read.table(infile, sep="\t", header=TRUE)
close(infile)

revCumRowSums <- function (x) {
    for (j in ncol(x):2) x[,j] <- rowSums(x[,1:j])
    return(x)
}

repro_self <- repro[,grepl("self_", names(repro))] |> revCumRowSums()
repro_ma <- repro[,grepl("ma_", names(repro))] |> revCumRowSums()
repro_pa <- repro[,grepl("pa_", names(repro))] |> revCumRowSums()
midp <- (repro_ma + repro_pa)/2
seg <- repro_self - midp

pop <- read.table(paste0(base, ".pop.tsv"), sep="\t", header=TRUE)
pop$age <- pop$age * params$DT
fix <- read.table(paste0(base, ".fix.tsv"), sep="\t", header=TRUE)

png(sprintf("%s_stratified.png", base), width=6.5, height=length(params$EPSILON)*1.5, pointsize=10, units='in', res=144)
{
    layout(matrix(seq_along(params$EPSILON), ncol=2, byrow=TRUE))
    par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
    for (j in seq_along(params$EPSILON)) {
        plot(midp[,j],
             seg[,j],
             pch=20, col=adjustcolor("black", 0.2),
             xlab=sprintf("midparent, effects <= %.3f", c(params$EPSILON, Inf)[j+1]),
             ylab="offspring - midparent"
        )
    }
}
dev.off()

plot_conditional_hists <- function (midp, seg) {
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
             xlim=c(-1,1) * lims)
        legend("bottomleft", title=sprintf("tail: %.2f%% < %.0f", 100 * mean(seg[ut] < -lims), -lims), legend='', bty='n')
        legend("bottomright", title=sprintf("tail: %.2f%% > %.0f", 100 * mean(seg[ut] > lims), lims), legend='', bty='n')
    }
    mtext("segregation noise", 1, line=2.5, cex=0.75)
}

for (j in seq_along(params$EPSILON)) {
    png(sprintf("%s_ind_%d.png", base, j), width=3.5, height=5, pointsize=10, units='in', res=144)
    plot_conditional_hists(midp[,j], seg[,j])
    dev.off()
}

plot_conditional_ratios <- function (midp, seg, do_legend=TRUE, title=NA, ...) {
    pbreaks <- pmax(0.01, pmin(0.99, seq(0, 1, length.out=51)))
    qbreaks <- c(-Inf, quantile(seg, pbreaks), Inf)
    total_xh <- hist(seg, breaks=qbreaks, plot=FALSE)
    uk <- seq(2, length(qbreaks)-1)
    par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
    plot(0, 0, type='n', xlim=range(qbreaks, finite=TRUE), ylim=c(0, 2), 
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

pdf(sprintf("%s_cond.pdf", base), width=9, height=6, pointsize=10)
par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
for (j in 1:ncol(seg)) {
    plot_conditional_ratios(midp[,j], seg[,j], title=sprintf("effects <= %.2f", c(params$EPSILON, Inf)[j+1]))
}
dev.off()

ft <- function (t, x) {
    out <- rep(NA, length(t))
    for (j in seq_along(t))  out[j] <- mean(cos(t[j] * x))
    return(out)
}

theory <- function (t, scale, gamma) {
    u <- gamma * t
    return( exp(2 * scale * ( 1 - u * pracma::Si(u) - cos(u) ) / (gamma*pi)) )
}

qtcauchy <- function(p, x, scale) {
    u <- pcauchy(-x, scale=scale);
    return( qcauchy(u + p * (1 - 2 * u), scale=scale) ) 
}

rtcauchy <- function (n, scale, gamma, M=1000) {
    out <- rep(NA, n)
    for (j in 1:n) {
        x <- qtcauchy(runif(M), gamma*M, scale=scale)
        out[j] <- mean(x)
    }
    return(out)
}

# check our simulation and theoretical FT agree
if (FALSE) {
    tt <- seq(0, 4, length.out=51)
    sim_scale <- 0.4; sim_gamma <- 0.1
    simvals <- rtcauchy(1e5, scale=sim_scale, gamma=sim_gamma)
    sim_ft <- ft(tt, simvals)
    layout(1:2)
    hist(simvals, breaks=40)
    plot(tt, sim_ft, type='l', xlab="t", ylab="fourier transform")
    lines(tt, theory(tt, scale=sim_scale, gamma=sim_gamma), lty=3, col='red')
    legend("topright", legend=c("simulation", "theory"), lty=c(1,3), col=c('black', 'red'))
}

tt <- seq(0, 4, length.out=51)
seg_var <- apply(seg, 2, var)
seg_ft <- do.call(cbind, lapply(1:ncol(seg), function (j) ft(tt, seg[,j]/sqrt(seg_var[j]))))

pdf(file=sprintf("%s_theory_observed_ft.pdf", base), width=9, height=6, pointsize=10)
layout(matrix(seq_along(params$EPSILON), nrow=2))
gvals <- seq(0.5, 4, length.out=11)
# variance is 2 * scale * gamma / pi
theory_ft <- do.call(cbind, lapply(gvals, function (g) { theory(tt, scale=g, gamma=pi / (2 * g)) }))
for (j in seq_along(params$EPSILON)) {
    plot(tt, seg_ft[,j], type='l', lwd=2,
            xlab='t', ylab='E[exp(itX)]')
    matlines(tt, theory_ft, lty=2, col=rainbow(10))
    lines(tt, exp(-tt^2/2), lty=3, col='black')
    legend("topright", lty=1:3, legend=c(paste('effects <=', c(params$EPSILON[-1], Inf)[j]), 'theory', 'gaussian'),
            col=c('black', 'red', 'black'))
    # legend("bottomleft", lty=2, col=rainbow(10)[seq_along(gvals)], legend=sprintf("g=%.3f", gvals))
}
dev.off()

# compare simulated distributions to observed
gvals <- seq(0.5, 4, length.out=11)
simvals <- do.call(cbind, lapply(gvals, function (g) sort(rtcauchy(1e5, scale=g, gamma=pi/(2*g)))))

png(file=sprintf("%s_sim_observed.png", base), width=9, height=6, pointsize=10, units='in', res=144)
par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
layout(matrix(seq_along(params$EPSILON), nrow=2))
for (j in seq_along(params$EPSILON)) {
    matplot(sort(sample(seg[,j], min(nrow(simvals), nrow(seg))))/sqrt(seg_var[j]), simvals, pch=20, xlab='observed', ylab='simulated')
    abline(0, 1)
}
dev.off()
