#!/usr/bin/Rscript

usage <- "
   make_figs.R (basename)
where basename.repro.tsv and basename.pop.tsv and basename.fix.tsv should be files,
and the 'repro' file should have a first line with parameters,
and remaining lines of the form
  (time) (trait values) (maternal trait values) (paternal trait values)
the 'pop' file should have 'trait' and 'age' columns,
and the 'fix' file should have 'time' and 'num_fixations' columns.
"

library(jsonlite)

args <- commandArgs(TRUE)
if (length(args) != 1) {
    stop(usage)
}

basename <- args[1]

infile <- file(paste0(basename, ".repro.tsv"), "r")
params <- fromJSON(readLines(infile, 1))
repro <- do.call(rbind, lapply(strsplit(readLines(infile), "\t"), function (x) {
           x <- as.numeric(x); n <- (length(x) - 1) / 3; data.frame(t=x[1], x=x[1+(1:n)], ma=x[1+n+(1:n)], pa=x[1+2*n+(1:n)])
}))
close(infile)

if (! "DT" %in% names(params)) params$DT <- 1.0

pop <- read.table(paste0(basename, ".pop.tsv"), sep="\t", header=TRUE)
pop$age <- pop$age * params$DT
fix <- read.table(paste0(basename, ".fix.tsv"), sep="\t", header=TRUE)

repro$midp <- (repro$ma + repro$pa)/2
repro$seg <- repro$x - repro$midp

# pdf(sprintf("%s_independence.pdf", basename), width=6.5, height=4, pointsize=10)
png(sprintf("%s_independence.png", basename), width=6.5, height=4, pointsize=10, units='in', res=144)
{
    layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), ncol=2), widths=c(1.5,3))
    par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
    # subset for plotting in a way that won't affect the plot
    # ut <- sample.int(nrow(repro), 3e4)
    # ut <- sort(unique(c(ut, which(
    #            pmax(rank(abs(repro$midp)), rank(abs(repro$seg))) > nrow(repro) * .99
    # ))))
    # # another way of doing this??
    # ut <- sample.int(nrow(repro),
    #             prob=pmax(rank(abs(repro$midp)), rank(abs(repro$seg)))^3,
    #             size=min(3e4, nrow(repro))
    # )
    ut <- TRUE
    for (q in c(.98, .998)) {
        lims <- quantile(abs(repro$midp), q)
        plot(
             repro$midp[ut],
             repro$seg[ut],
             pch=20, cex=0.5, col=adjustcolor('black', 0.2), asp=1,
             xlab='midparent trait',
             ylab='segregation noise',
             xlim=c(-1,1) * lims,
             ylim=c(-1,1) * lims
        )
    }
    lims <- quantile(abs(repro$seg), 0.98)
    for (q in c(.1, .5, .9)) {
        ut <- (abs(rank(repro$midp)/nrow(repro) - q) < 0.025)
        outl <- (repro$seg > -lims) & (repro$seg < lims)
        hist(repro$seg[ut & outl], breaks=100, main=sprintf("midparent %0.0f%%", 100*q),
             # xlab="segregation noise",
             xlab='',
             xlim=c(-1,1) * lims)
        legend("bottomleft", title=sprintf("tail: %.2f%% < %.0f", 100 * mean(repro$seg[ut] < -lims), -lims), legend='', bty='n')
        legend("bottomright", title=sprintf("tail: %.2f%% > %.0f", 100 * mean(repro$seg[ut] > lims), lims), legend='', bty='n')
    }
    mtext("segregation noise", 1, line=2.5, cex=0.75)
}
dev.off()

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

pdf(sprintf("%s_dist.pdf", basename), width=6, height=2, pointsize=10)
layout(t(1:2))
par(mar=c(4,4,1,1), mgp=c(2.5,1,0))
    plot_tail(repro$midp, xlab='midparent')
    plot_tail(repro$seg, xlab='segregation noise')
dev.off()
