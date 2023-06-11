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

phenotypes <- read.table("phenotypes.both_sexes.tsv", stringsAsFactors=FALSE, sep="\t", quote="", header=TRUE, comment.char="")

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
results <- read.table("snp_results_10.tsv", stringsAsFactors=FALSE)
results <- subset(results, !is.na(slope))
results$name <- manifest$Phenotype.Description[match(results$code, manifest$Phenotype.Code)]
results$name[results$name == "NA"] <- NA
results <- subset(results, !grepl("Job code", results$name))
results <- subset(results, !grepl("Type of", results$name)) # seems to be all diet-related
results <- subset(results, !grepl("Own or rent", results$name))
results <- subset(results, !grepl("employment", results$name))
results <- subset(results, !grepl("Milk type", results$name))
results <- subset(results, !grepl("Spread type", results$name))
results <- subset(results, !grepl("Bread type", results$name))
results <- subset(results, !grepl("Types of transport", results$name))

pheno_names <- c("phenotype", "n_non_missing", "n_missing", "n_controls", "n_cases")
pheno <- phenotypes[,pheno_names]
names(pheno) <- c("code", pheno_names[-1])
stopifnot(all(results$code %in% pheno$code))
results <- merge(results, pheno, all.x=TRUE)
# very common things (with > 10000 cases) have Something Different going on
results$pruned <- !(!is.na(results$n_cases) & results$n_non_missing > 300000 & results$n_cases < 10000)

categories <- c(
    "Non-cancer illness code",
    "neoplasm",
    "ICD10",
    "other"
)
results$category <- ifelse(
    grepl(categories[1], results$name), categories[1],
    ifelse(
        grepl(categories[2], results$name), categories[2],
        ifelse(
            grepl(categories[3], results$name), categories[3], categories[4]
        )
    )
)

for (cname in categories) {
    pdf(sprintf("%s_results_10.pdf", gsub(" ", "_", cname)), width=7, height=6)
    with(subset(results, category == cname), {
        layout(matrix(c(1,2,1,3,1,4), nrow=2))
        hist(slope, breaks=30, xlab='tail exponent', main=cname)
        plot(nsnps, slope, pch=20, cex=0.5, xlab='number of snps', ylab='tail exponent',
             col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
        plot(n_cases, slope, pch=20, cex=0.5, xlab='number of cases', ylab='tail exponent',
             col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
        plot(n_non_missing, slope, pch=20, cex=0.5, xlab='number of non-missing subjects', ylab='tail exponent',
             col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
    })
    dev.off()
}

pdf("unfiltered_results_10.pdf", width=7, height=6, pointsize=10)
with(results, {
    layout(matrix(c(1,2,1,3,1,4), nrow=2))
    hist(slope, breaks=30, xlab='tail exponent', main='')
    plot(nsnps, slope, pch=20, cex=0.5, xlab='number of snps', ylab='tail exponent',
         col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
    plot(n_cases, slope, pch=20, cex=0.5, xlab='number of cases', ylab='tail exponent',
         col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
    plot(n_non_missing, slope, pch=20, cex=0.5, xlab='number of non-missing subjects', ylab='tail exponent',
         col=adjustcolor(ifelse(pruned, "red", "black"), 0.5))
})
dev.off()

pdf(file="results_10.pdf", width=6, height=4, pointsize=10)
layout(matrix(c(1,2,1,3), nrow=2))
with(subset(results, !pruned), {
    par(mar=c(4, 4, 1, 1), mgp=c(2.1,1,0))
    hist(slope, breaks=50, xlab='tail exponent', main='')
mtext("(A)", 3, adj=-0.1, line=0)
    plot(nsnps, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5), log='x', xlab="number of SNPs", ylab='tail exponent')
mtext("(B)", 3, adj=-0.25, line=0)
    plot(n_cases, slope, pch=20, cex=0.5, col=adjustcolor('black', 0.5), log='x', xlab="number of cases", ylab='tail exponent')
mtext("(C)", 3, adj=-0.25, line=0)
})
dev.off()

pdf(file="slopes_by_category.pdf", width=6, height=9, pointsize=10)
layout(seq_along(categories))
breaks <- hist(results$slope[!results$pruned], breaks=30, plot=FALSE)$breaks
for (cname in categories) {
    with(subset(results, !pruned & category == cname), {
        hist(slope, breaks=breaks, main=cname)
    })
}
dev.off()

clean_name <- function (x) {
    clean <- c("Non-cancer illness code, self-reported: ", 
               "Diagnoses - main ICD10: ",
               "Cancer code, self-reported: "
    )
    for (u in clean) {
        x <- gsub(u, "", x)
    }
    return(x)
}

plot_code <- function (code, lab) {
    ij <- match(code, results$code)
    if (is.na(ij)) stop(sprintf("Code %s not in results.", code))
    code_dir <- file.path("phenotypes", code)
    xlm_file <- file.path(code_dir, paste0("tail_10_", "xlm.RData"))
    if (!file.exists(xlm_file)) stop(sprintf("File %s does not exist.", xlm_file))
    load(xlm_file)
    logx <- xlm$model[["log(x)"]]
    logy <- xlm$model[["log(tail_prob)"]]
    plot(logx, logy, ylab='log(prob > t)', xlab='log(t)', type='l')#, sub=code)
    cnp <- strsplit(clean_name(results$name[ij]), " ")[[1]]
    title(paste(c(lab, cnp[1:min(5, length(cnp))]), collapse=" "), cex.main=0.9)
    if (length(cnp) > 5)
        title(paste(cnp[6:length(cnp)], collapse=" "), line=0.2, cex.main=0.9)
    abline(coef(xlm), col='red')
    legend("topright", lty=c(1, NA, NA), cex=0.8, col=c('red', NA, NA),
           legend=c(sprintf("y = %0.2f %0.2f x", coef(xlm)[1], coef(xlm)[2]),
                    sprintf("num snps = %d", results$nsnps[ij]),
                    sprintf("num cases = %d", results$n_cases[ij])
           )
    )
}

set.seed(8)
examples <- lapply(categories, function (cname) {
       with(subset(results, category==cname & !pruned & !is.na(name)), {
            xc <- cut(slope, c(-100,-2,-1.5,0))
            tapply(code, xc, sample, 1)
       })
})

pdf(file="examples.pdf", width=6.5, height=8, pointsize=10)
layout(matrix(1:12, nrow=4, byrow=TRUE))
par(mar=c(4, 4, 2, 1))
for (j in 1:12) {
    cname <- unlist(examples)[j]
    plot_code(cname, paste0("(", letters[j], ")"))
}
dev.off()
