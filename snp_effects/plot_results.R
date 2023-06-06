#!/usr/bin/env Rscript

## First run with prop 25%
results <- read.table("snp_results_25.tsv")
pdf(file="results_25.pdf", width=6, height=8, pointsize=10)
layout((1:2))
with(subset(results, slope > -4), {
    hist(slope, breaks=50, xlab='exponent', main="tail proportion: 25%")
    plot(nsnps, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5))
})
dev.off()


## Next run with prop 10%
results <- read.table("snp_results_10.tsv")
pdf(file="results_10.pdf", width=6, height=8, pointsize=10)
layout((1:2))
with(subset(results, slope > -4), {
    hist(slope, breaks=50, xlab='exponent', main="tail proportion: 10%")
    plot(nsnps, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5))
})
dev.off()

