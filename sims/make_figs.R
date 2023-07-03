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
source("helpers.R")

args <- commandArgs(TRUE)
if (length(args) != 1) {
    stop(usage)
}

basename <- args[1]

infile <- file(paste0(basename, ".repro.tsv"), "r")
params <- fromJSON(readLines(infile, 1))
repro <- read.table(infile, sep="\t", header=TRUE)
close(infile)

repro_self <- repro[,grepl("self_", names(repro))] |> revCumRowSums()
repro_ma <- repro[,grepl("ma_", names(repro))] |> revCumRowSums()
repro_pa <- repro[,grepl("pa_", names(repro))] |> revCumRowSums()
midp <- (repro_ma + repro_pa)/2
seg <- repro_self - midp

if (! "DT" %in% names(params)) params$DT <- 1.0

pop <- read.table(paste0(basename, ".pop.tsv"), sep="\t", header=TRUE)
pop$age <- pop$age * params$DT
fix <- read.table(paste0(basename, ".fix.tsv"), sep="\t", header=TRUE)

# repro <- do.call(rbind, lapply(strsplit(readLines(infile), "\t"), function (x) {
#            x <- as.numeric(x); n <- (length(x) - 1) / 3; data.frame(t=x[1], x=x[1+(1:n)], ma=x[1+n+(1:n)], pa=x[1+2*n+(1:n)])
# }))
# repro$midp <- (repro$ma + repro$pa)/2
# repro$seg <- repro$x - repro$midp

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
    j <- ncol(midp)
    for (q in c(.98, .998)) {
        lims <- quantile(abs(midp[,j]), q)
        plot(
             midp[ut,j],
             seg[ut,j],
             pch=20, cex=0.5, col=adjustcolor('black', 0.2), asp=1,
             xlab='midparent trait',
             ylab='segregation noise',
             xlim=c(-1,1) * lims,
             ylim=c(-1,1) * lims
        )
    }
    lims <- quantile(abs(seg[,j]), 0.98)
    for (q in c(.1, .5, .9)) {
        ut <- (abs(rank(midp[,j])/nrow(midp) - q) < 0.025)
        outl <- (seg[,j] > -lims) & (seg[,j] < lims)
        hist(seg[ut & outl,j], breaks=100, main=sprintf("midparent %0.0f%%", 100*q),
             # xlab="segregation noise",
             xlab='',
             xlim=c(-1,1) * lims)
        legend("bottomleft", title=sprintf("tail: %.2f%% < %.0f", 100 * mean(seg[ut,j] < -lims), -lims), legend='', bty='n')
        legend("bottomright", title=sprintf("tail: %.2f%% > %.0f", 100 * mean(seg[ut,j] > lims), lims), legend='', bty='n')
    }
    mtext("segregation noise", 1, line=2.5, cex=0.75)
}
dev.off()

png(sprintf("%s_dist.png", basename), width=6, height=2, pointsize=10, units='in', res=144)
layout(t(1:2))
par(mar=c(4,4,1,1), mgp=c(2.5,1,0))
    j <- ncol(midp)
    plot_tail(midp[,j], xlab='midparent')
    plot_tail(seg[,j], xlab='segregation noise')
dev.off()
