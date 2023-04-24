#!/usr/bin/Rscript


usage <- "
   parse_output.R (basename) (number of time bins) 
where basename.repro.tsv and bsaename.pop.tsv and basename.fix.tsv should be files,
and the 'repro' file should have a first line with parameters,
and remaining lines of the form
  (time) (trait values) (maternal trait values) (paternal trait values)
the 'pop' file should have 'trait' and 'age' columns,
and the 'fix' file should have 'time' and 'num_fixations' columns.
"

library(jsonlite)

args <- commandArgs(TRUE)
if (length(args) != 2) {
    stop(usage)
}

basename <- args[1]
num_plots <- as.numeric(args[2])

infile <- file(paste0(basename, ".repro.tsv"), "r")
params <- fromJSON(readLines(infile, 1))
repro <- do.call(rbind, lapply(strsplit(readLines(infile), "\t"), function (x) {
           x <- as.numeric(x); n <- (length(x) - 1) / 3; data.frame(t=x[1], x=x[1+(1:n)], ma=x[1+n+(1:n)], pa=x[1+2*n+(1:n)])
}))
close(infile)

pop <- read.table(paste0(basename, ".pop.tsv"), sep="\t", header=TRUE)
pop$age <- pop$age * params$DT
fix <- read.table(paste0(basename, ".fix.tsv"), sep="\t", header=TRUE)

pdf(sprintf("%s.fixrate.pdf", basename), width=6, height=6, pointsize=10)
layout(1:2)
    plot(fix$time, fix$num_fixations, type='l', xlab='time', ylab='number of fixations')
    hist(diff(fix$num_fixations), xlab='number of fixations/tick', breaks=40)
dev.off()

pdf(sprintf("%s.ages.pdf", basename), width=6, height=8, pointsize=10)
layout(1:3)
    hist(pop$age, main=toJSON(params), xlab="age", breaks=30)
    hist(pop$trait, xlab='trait', main='trait', breaks=30)
    hist(log(abs(pop$trait)), xlab='log|trait|', main='log|trait|', breaks=30)
dev.off()

repro$seg <- repro$x - (repro$ma + repro$pa)/2
repro$plot <- cut(repro$t, num_plots)

if (params$type == "normal") {
    qfun <- qnorm
} else if (params$type == "cauchy") {
    qfun <- qcauchy
}

pdf(sprintf("%s.seg.pdf", basename), width=8, height=2*num_plots, pointsize=10)
layout(matrix(1:(3*num_plots), ncol=3, byrow=TRUE))
xlim <- range(repro$seg)
xh <- hist(repro$seg, breaks=200, plot=FALSE)
for (k in levels(repro$plot)) {
    title <- if (k == levels(repro$plot)[1]) toJSON(params) else ""
    sub <- repro[repro$plot == k,]
    hist(sub$seg, breaks=xh$breaks,
         main=sprintf("times %0.2f-%0.2f", min(sub$t), max(sub$t)),
         xlab="trait - midparent",
         xlim=xlim
    )
    qqplot(qfun(ppoints(sub$seg)), sub$seg, ylab="trait value", xlab="cauchy quantiles", main=title)
    qqline(sub$seg, distribution = function(p) qfun(p), probs = c(0.1, 0.6), col = 2)
    x <- sort(abs(sub$seg))
    tail_prob <- 1 - seq_along(x)/(length(x) + 1)
    k <- max(sum(x==0) + 1, floor(length(x) * 3/4)) : length(x)
    plot(log(x[k]), log(tail_prob[k]), xlab="log trait value", ylab="log P(X>x)")
    ab <- coef(lm(log(tail_prob[k]) ~ log(x[k])))
    abline(ab, col='red')
    legend("topright", lty=1, legend=sprintf("y = %0.2f x + %0.2f", ab[2], ab[1]))
}
dev.off()
