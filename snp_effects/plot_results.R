#!/usr/bin/env Rscript

manifest <- subset(read.csv("gwas_manifest.csv", header=TRUE, na.string="N/A", stringsAsFactors=FALSE),
                   Sex=="both_sexes" & !is.na(Phenotype.Code) & !grepl("_irnt", Phenotype.Code)
                   & !grepl("^20003_", Phenotype.Code) # drug treatments
                   & !grepl("22617_", Phenotype.Code) # job codes
                   & !grepl("22601_", Phenotype.Code) # also job codes
            )
# remove doubles
doubles <- names(table(manifest$Phenotype.Code))[table(manifest$Phenotype.Code) == 2]
manifest <- manifest[-(match(doubles,manifest$Phenotype.Code)),]

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
results$name <- manifest$Phenotype.Description[match(results$code, manifest$Phenotype.Code)]
results <- subset(results, !grepl("Job code", results$name))
results <- subset(results, !grepl("Type of", results$name)) # seems to be all diet-related

pdf(file="results_10.pdf", width=6, height=8, pointsize=10)
layout((1:2))
with(subset(results, slope > -4), {
    hist(slope, breaks=50, xlab='exponent', main="tail proportion: 10%")
    plot(nsnps, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5))
})
dev.off()

categories <- c(
    "Non-cancer illness code",
    "Cancer code",
    "neoplasm",
    "ICD10"
)

pdf("cat_results_10.pdf", width=6, height=9)
for (cname in categories) {
    with(subset(results, grepl(cname, name)), {
        layout((1:3))
        hist(slope, breaks=50, xlab='exponent', main=cname)
        hist(slope[slope > -4], breaks=50, xlab='exponent', main="tail proportion: 10%")
        plot(nsnps, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5))
    })
}
dev.off()
