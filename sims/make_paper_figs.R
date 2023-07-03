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

# 1. Plots of median trait value over time and final trait distribution
#  for neutral Normal and Cauchy examples
fname <- "neutral_trait_traces.pdf"
examples <- lapply(list(
                 "neutral Normal" = "sim_neutral_normal_831",
                 "neutral Cauchy" = "sim_neutral_cauchy_831"
    ), load_files)
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

# 2. Plots of median trait value over time and final trait distribution
#  for selected Normal and Cauchy examples
fname <- "selected_trait_traces.pdf"
examples <- lapply(list(
                 "selected Normal" = "sim_normal_707",
                 "selected Cauchy" = "sim_cauchy_707"
    ), load_files)
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


