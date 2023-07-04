#!/usr/bin/Rscript
library(jsonlite)
source("helpers.R")

load_files <- function (base) {
    infile <- file(paste0(base, ".repro.tsv"), "r")
    params <- fromJSON(readLines(infile, 1))
    repro <- read.table(infile, sep="\t", header=TRUE)
    close(infile)
    if (! "DT" %in% names(params)) params$DT <- 1.0
    pop <- read.table(paste0(base, ".pop.tsv"), sep="\t", header=TRUE)
    pop$age <- pop$age * params$DT
    fix <- read.table(paste0(base, ".fix.tsv"), sep="\t", header=TRUE)
    repro_self <- repro[,grepl("self_", names(repro))] |> revCumRowSums()
    repro_ma <- repro[,grepl("ma_", names(repro))] |> revCumRowSums()
    repro_pa <- repro[,grepl("pa_", names(repro))] |> revCumRowSums()
    midp <- (repro_ma + repro_pa)/2
    seg <- repro_self - midp
    return(list(
                basename=base,
                params=params,
                repro_self=repro_self,
                midp=midp,
                seg=seg,
                pop=pop,
                fix=fix
    ))
}

plot_seg_noise <- function (examples, labels, seg_col, main_append="", ...) {
    layout(matrix(1:6, nrow=3), heights=c(1.5,1, 1))
    for (x in names(examples)) {
        j <- if (!missing(seg_col)) seg_col else ncol(examples[[x]][['seg']])
        main <- paste0(x, 
                       if (missing(seg_col)) ""
                       else
                       sprintf(", effects <= %.2f", c(examples[[x]][['params']][['EPSILON']],Inf)[j+1])
        )
        this_seg <- examples[[x]][['seg']][,j]
        this_midp <- examples[[x]][['midp']][,j]
        par(mar=c(4,4,2,0)+.1, mgp=c(2.5,1,0))
            plot(this_midp, this_seg, pch=20, col=adjustcolor("black", 0.5), asp=1,
                 xlab="midparent trait", ylab="segregation noise",
                 main=main)
        mtext(labels[[x]][1], 3, adj=-0.0, line=0.2)
        par(mar=c(4,4,0,0)+.1, mgp=c(2.5,1,0))
            hist(this_seg, breaks=40, main='',
                 xlab='segreation noise')
        mtext(labels[[x]][2], 3, adj=-0.0, line=0.2)
            plot_conditional_ratios(this_midp, this_seg, do_legend=(x == names(examples)[1]), ...)
        mtext(labels[[x]][3], 3, adj=-0.0, line=0.2)
    }
}


###########
# Neutral examples
examples <- lapply(list(
                 "neutral Normal" = "sim_neutral_normal_831",
                 "neutral Cauchy" = "sim_neutral_cauchy_831"
    ), load_files)

# 1. Plots of median trait value over time and final trait distribution
fname <- "neutral_trait_traces.pdf"
labels <- list(c("(A)", "(B)"), c("(C)", "(D)"))
names(labels) <- names(examples)

pdf(file=fname, width=6.5, height=3, pointsize=10)
layout(matrix(1:4, nrow=2), heights=c(1.5,1))
for (x in names(examples)) {
    par(mar=c(4,4,2,0)+.1, mgp=c(2.5,1,0))
        plot(med_trait ~ time, data=examples[[x]][['fix']], type='l',
             main=x, xlab='time', ylab='median trait value')
    mtext(labels[[x]][1], 3, adj=-0.0, line=0.2)
    par(mar=c(4,4,0,0)+.1, mgp=c(2.5,1,0))
        hist(examples[[x]][['pop']]$trait, breaks=100,
             xlab='trait value', main='')
    mtext(labels[[x]][2], 3, adj=-0.0, line=0.2)
}
dev.off()

# 3. Segregation noise:
# (distrn of seg noise), (seg noise vs midparent), and (Radon-Nicodym deriv of conditional seg noise for a few values of midparent)

fname <- "neutral_seg_noise.png"
labels <- list(c("(A)", "(B)", "(C)"), c("(D)", "(E)", "(F)"))
names(labels) <- names(examples)

png(file=fname, width=6.5, height=6, pointsize=10, units='in', res=288)
plot_seg_noise(examples, labels)
dev.off()

# 4. Segregation noise for small effects
# (distrn of seg noise), (seg noise vs midparent), and (Radon-Nicodym deriv of conditional seg noise for a few values of midparent)
# ONLY for small effects

fname <- "neutral_seg_noise_small.png"
labels <- list(c("(A)", "(B)", "(C)"), c("(D)", "(E)", "(F)"))
names(labels) <- names(examples)

png(file=fname, width=6.5, height=6, pointsize=10, units='in', res=288)
plot_seg_noise(examples, labels, seg_col=1)
dev.off()


##########
#  Examples with selection
examples <- lapply(list(
                 "selected Normal" = "sim_normal_707",
                 "selected Cauchy" = "sim_cauchy_707"
    ), load_files)

# 1. Plots of median trait value over time and final trait distribution
fname <- "selected_trait_traces.pdf"
labels <- list(c("(A)", "(B)"), c("(C)", "(D)"))
names(labels) <- names(examples)

pdf(file=fname, width=6.5, height=3, pointsize=10)
layout(matrix(1:4, nrow=2), heights=c(1.5,1))
for (x in names(examples)) {
    par(mar=c(4,4,2,0)+.1, mgp=c(2.5,1,0))
        plot(med_trait ~ time, data=examples[[x]][['fix']], type='l',
             main=x, xlab='time', ylab='median trait value')
    mtext(labels[[x]][1], 3, adj=-0.0, line=0.2)
    par(mar=c(4,4,0,0)+.1, mgp=c(2.5,1,0))
        hist(examples[[x]][['pop']]$trait, breaks=100,
             xlab='trait value', main='')
    mtext(labels[[x]][2], 3, adj=-0.0, line=0.2)
}
dev.off()

# Selected

# 3. Segregation noise:
# (distrn of seg noise), (seg noise vs midparent), and (Radon-Nicodym deriv of conditional seg noise for a few values of midparent)

fname <- "selected_seg_noise.png"
labels <- list(c("(A)", "(B)", "(C)"), c("(D)", "(E)", "(F)"))
names(labels) <- names(examples)

png(file=fname, width=6.5, height=6, pointsize=10, units='in', res=288)
plot_seg_noise(examples, labels, w=0.1)
dev.off()


# 4. Segregation noise for small effects
# (distrn of seg noise), (seg noise vs midparent), and (Radon-Nicodym deriv of conditional seg noise for a few values of midparent)
# ONLY for small effects

fname <- "selected_seg_noise_small.png"
labels <- list(c("(A)", "(B)", "(C)"), c("(D)", "(E)", "(F)"))
names(labels) <- names(examples)

png(file=fname, width=6.5, height=6, pointsize=10, units='in', res=288)
plot_seg_noise(examples, labels, seg_col=1)
dev.off()

