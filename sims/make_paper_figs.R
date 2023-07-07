#!/usr/bin/Rscript
library(jsonlite)
source("helpers.R")

load_params <- function (base) {
    infile <- file(paste0(base, ".repro.tsv"), "r")
    params <- fromJSON(readLines(infile, 1))
    close(infile)
    if (! "DT" %in% names(params)) params$DT <- 1.0
    return(params)
}

load_fix <- function (base) {
    fix <- read.table(paste0(base, ".fix.tsv"), sep="\t", header=TRUE)
    return(fix)
}

load_pop <- function (base) {
    params <- load_params(base)
    pop <- read.table(paste0(base, ".pop.tsv"), sep="\t", header=TRUE)
    pop$age <- pop$age * params$DT
    return(pop)
}

load_seg_midp <- function (base, max_n=as.integer(1e6)) {
    infile <- file(paste0(base, ".repro.tsv"), "r")
    params <- fromJSON(readLines(infile, 1))
    repro <- read.table(infile, sep="\t", header=TRUE, nrows=max_n)
    close(infile)
    if (! "DT" %in% names(params)) params$DT <- 1.0
    repro_self <- repro[,grepl("self_", names(repro))] |> revCumRowSums()
    repro_ma <- repro[,grepl("ma_", names(repro))] |> revCumRowSums()
    repro_pa <- repro[,grepl("pa_", names(repro))] |> revCumRowSums()
    midp <- (repro_ma + repro_pa)/2
    seg <- repro_self - midp
    return(list(midp=midp, seg=seg))
}

plot_trait_traces <- function (examples, labels) {
    layout(matrix(1:(2*length(examples)), nrow=2), heights=c(1.5,1))
    for (x in names(examples)) {
        fix <- load_fix(examples[[x]])
        pop <- load_pop(examples[[x]])
        par(mar=c(4,4,2,0)+.1, mgp=c(2.5,1,0))
            plot(med_trait ~ time, data=fix, type='l',
                 main=x, xlab='time', ylab='median trait value')
        mtext(labels[[x]][1], 3, adj=-0.0, line=0.2)
        par(mar=c(4,4,0,0)+.1, mgp=c(2.5,1,0))
            hist(pop$trait, breaks=100,
                 xlab='trait value', main='')
        mtext(labels[[x]][2], 3, adj=-0.0, line=0.2)
        rm(fix)
        rm(pop)
    }
}

plot_seg_noise <- function (examples, labels, seg_col, main_append="", ...) {
    layout(matrix(1:(3*length(examples)), nrow=3), heights=c(1.5,1, 1))
    for (x in names(examples)) {
        params <- load_params(examples[[x]])
        sm <- load_seg_midp(examples[[x]])
        j <- if (!missing(seg_col)) seg_col else ncol(sm[['seg']])
        main <- paste0(x, 
                       if (missing(seg_col)) ""
                       else
                       sprintf(", effects <= %.2f", c(params[['EPSILON']],Inf)[j+1])
        )
        this_seg <- sm[['seg']][,j]
        this_midp <- sm[['midp']][,j]
        par(mar=c(4,4,2,0)+.1, mgp=c(2.5,1,0))
            plot(this_midp, this_seg, pch=20, col=adjustcolor("black", 0.5), asp=1,
                 xlab="midparent trait", ylab="segregation noise",
                 main=main)
        mtext(labels[[x]][1], 3, adj=-0.0, line=0.2)
        par(mar=c(4,4,0,0)+.1, mgp=c(2.5,1,0))
            hist(this_seg, breaks=40, main='',
                 xlab='segregation noise')
        mtext(labels[[x]][2], 3, adj=-0.0, line=0.2)
            plot_conditional_ratios(this_midp, this_seg, do_legend=(x == names(examples)[1]), ...)
        mtext(labels[[x]][3], 3, adj=-0.0, line=0.2)
    }
}


###########
# Neutral examples
examples <- list(
                 # "neutral Normal" = "sim_neutral_normal_831",
                 # "neutral Cauchy" = "sim_neutral_cauchy_831"
                 "neutral Normal" = "sim_neutral_normal_8135",
                 "neutral Cauchy" = "sim_neutral_cauchy_8135"
)

# 1. Plots of median trait value over time and final trait distribution
fname <- "neutral_trait_traces.pdf"
labels <- list(c("(A)", "(B)"), c("(C)", "(D)"))
names(labels) <- names(examples)

pdf(file=fname, width=6.5, height=3, pointsize=10)
plot_trait_traces(examples, labels)
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
examples <- list(
                 "selected Normal" = "sim_normal_707",
                 "selected Cauchy" = "sim_cauchy_707"
)

# 1. Plots of median trait value over time and final trait distribution
fname <- "selected_trait_traces.pdf"
labels <- list(c("(A)", "(B)"), c("(C)", "(D)"))
names(labels) <- names(examples)

pdf(file=fname, width=6.5, height=3, pointsize=10)
plot_trait_traces(examples, labels)
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


# 5. with stable(3/2)
examples <- list(
                 "neutral Normal" = "sim_neutral_normal_8135",
                 "neutral Stable(3/2)" = "sim_neutral_stable_8135",
                 "neutral Cauchy" = "sim_neutral_cauchy_8135"
)

fname <- "neutral_trait_traces2.pdf"
labels <- list(c("(A)", "(B)"), c("(C)", "(D)"), c("(E)", "(F)"))
names(labels) <- names(examples)

pdf(file=fname, width=6.5, height=3, pointsize=10)
plot_trait_traces(examples, labels)
dev.off()

# 6. Segregation noise:

fname <- "neutral_seg_noise2.png"
labels <- list(c("(A)", "(B)", "(C)"), c("(D)", "(E)", "(F)"), c("(G)", "(H)", "(I)"))
names(labels) <- names(examples)

png(file=fname, width=8.5, height=6, pointsize=10, units='in', res=288)
plot_seg_noise(examples, labels)
dev.off()

