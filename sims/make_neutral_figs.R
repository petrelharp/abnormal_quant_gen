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

if (!exists("basename")) {
    args <- commandArgs(TRUE)
    if (length(args) != 1) {
        stop(usage)
    }

    basename <- args[1]
}

infile <- file(paste0(basename, ".repro.tsv"), "r")
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

pop <- read.table(paste0(basename, ".pop.tsv"), sep="\t", header=TRUE)
pop$age <- pop$age * params$DT
fix <- read.table(paste0(basename, ".fix.tsv"), sep="\t", header=TRUE)

png(sprintf("%s_stratified.png", basename), width=6.5, height=length(params$EPSILON)*1.5, pointsize=10, units='in', res=144)
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

plot_independence <- function (midp, seg) {
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
    png(sprintf("%s_ind_%d.png", basename, j), width=3.5, height=5, pointsize=10, units='in', res=144)
    plot_independence(midp[,j], seg[,j])
    dev.off()
}

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

# compare to Gaussian
png(file=sprintf("%s_theory_observed.png", basename), width=6.5, height=8, pointsize=10, units='in', res=144)
matplot(tt, seg_ft, type='l', lty=1, col=rainbow(10)[1:ncol(seg_ft)],
        xlab='t', ylab='E[exp(itX)]')
lines(tt, exp(-tt^2/2), lty=2)
legend("topright", lty=c(rep(1, ncol(seg_ft)), 2), col=c(rainbow(10)[1:ncol(seg_ft)], 'black'),
       legend=c(sprintf("effects <= %.3f", c(params$EPSILON[-1], Inf)), 'Gaussian'))
dev.off()

# which_j <- 3
# png(file=sprintf("%s_theory_observed.png", basename), width=6.5, height=8, pointsize=10, units='in', res=144)
# layout(1:2)
# svals <- exp(seq(log(0.01), log(10), length.out=11))
#     par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
# plot(tt, seg_ft, type='l', xlab="t", ylab="fourier transform")
#     for (g in gvals) {
#         lines(tt, theory(tt, var(seg[,which_j]), g), lty=3, col=rainbow(20)[match(g, gvals)])
#     }
#     legend("topright", lty=c(1, rep(3, length(gvals))), lwd=2, col=c('black', rainbow(20)[seq_along(gvals)]), legend=c('observed', sprintf("gamma = %.3f", gvals)))
# 
# ghat <- 0.013
# shat <- sqrt(var(seg[,which_j]) * pi / (ghat * 2))
# simvals <- rtcauchy(1e5, scale=shat, gamma=ghat)
# plot(sort(sample(seg[,which_j], length(simvals))), sort(simvals), xlab='observed', ylab='simulated')
# abline(0, 1)
# dev.off()
