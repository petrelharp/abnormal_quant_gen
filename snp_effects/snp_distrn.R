library(data.table)

# The distribution we are drawing X from is:
#  1. pick a SNP at random
#  2. flip a coin; with prob 1-minor_AF return 0
#  3. otherwise return beta
# and so with N SNPs
#   P(X > x)  = sum(minor_AF[beta > x])/N
# We are looking at abs(beta).

pval_cutoff <- 1e-8
min_num <- 20
# max_prop <- 0.25 ## FIRST VERSION WAS THIS
max_prop <- 0.1

filebase <- sprintf("tail_%d_", floor(100*max_prop))

estim_alpha <- function (gwas, num_breaks=1001) {
    # For a sequence of value xi,
    # fits the linear model
    #  log(pxi) ~ log(xi)
    # where
    #  pxi = prob( a randomly chosen locus has effect above xi )
    # to only the 'tail', defined as those xi such that there
    # are at least min_num loci with effect above xi
    # and at most a proportion max_prop of loci with effect above xi.
    xx <- c(seq(0, 2 * quantile(abs(gwas$beta), 0.99), length.out=num_breaks), Inf)
    disc_beta <- cut(abs(gwas$beta), xx)
    num_points <- table(disc_beta)
    dpx <- tapply(gwas$minor_AF, disc_beta, sum)
    dpx[is.na(dpx)] <- 0
    p_zero <- mean(gwas$minor_AF)
    tail_prob <- p_zero - cumsum(dpx) / nrow(gwas)
    tail_nums <- nrow(gwas) - cumsum(num_points)
    x <- xx[-1]
    use_these <- (is.finite(log(tail_prob))
                      & tail_nums > min_num
                      & tail_nums < max_prop * nrow(gwas))
    if (sum(use_these) > min_num) {
        xlm <- lm(log(tail_prob) ~ log(x), subset=use_these)
    } else {
        xlm <- NULL
    }
    return(xlm)
}

manifest <- subset(read.csv("gwas_manifest.csv", header=TRUE, na.string="N/A", stringsAsFactors=FALSE),
                   Sex=="both_sexes" & !is.na(Phenotype.Code) & !grepl("_irnt", Phenotype.Code)
                   & !grepl("^20003_", Phenotype.Code) # drug treatments
                   & !grepl("22617_", Phenotype.Code) # job codes
                   & !grepl("22601_", Phenotype.Code) # also job codes
            )
# remove doubles
doubles <- names(table(manifest$Phenotype.Code))[table(manifest$Phenotype.Code) == 2]
manifest <- manifest[-(match(doubles,manifest$Phenotype.Code)),]
translation <- read.csv("gwas_manifest_lookup.csv", header=FALSE)

# gwas_file <- "100002_raw.gwas.imputed_v3.female.tsv.bgz"
# gwas_file <- "3062_irnt.gwas.imputed_v3.both_sexes.tsv.bgz"

for (k in 1:nrow(manifest)) {
    code <- manifest$Phenotype.Code[k]
    if (!dir.exists(code)) dir.create(code)
    basename <- file.path(code, filebase)
    donefile <- paste0(basename, "done")
    xlm_file <- paste0(basename, "xlm.RData")
    plot_file <- paste0(basename, "plot.pdf")

    gwas_file <- file.path(code, manifest$File[k])
    if (!file.exists(gwas_file)) {
        wgot <- system(paste("cd", code, "&&", manifest$wget.command[k]))
        if (wgot != 0) {
            cat("Failed to get", code, "\n")
            next
        }
    }

    if (file.exists(donefile)) {
        cat("Already did", code, "\n")
        next
    }
    header <- scan(gwas_file, nlines=1, what='')
    pval_column <- which(header == "pval")
    read_cmd <- sprintf("zcat %s | head -n 1; zcat %s | awk '$%d < %0.16f'", gwas_file, gwas_file, pval_column, pval_cutoff)
    gwas <- fread(cmd=read_cmd)
    stopifnot(all(gwas$pval < pval_cutoff))
    if (ncol(gwas) == 0) {
        cat("Failed to load", code, "\n")
        next
    }
    gwas <- subset(gwas, !is.na(beta) & !is.na(minor_AF) & (beta != 0) & (minor_AF != 0))
    if (nrow(gwas) <= min_num) {
        cat("Not enough hits in", code, "\n")
        next
    }

    xlm <- estim_alpha(gwas)
    if (is.null(xlm)) {
        cat("Failed to do regression for", code, "\n")
        next
    }
    save(xlm, file=xlm_file)

    logy <- xlm$model[["log(tail_prob)"]]
    logx <- xlm$model[["log(x)"]]
    plotmain <- code
    if (code %in% translation[,1]) {
        plotmain <- translation[match(code, translation[,1]), 2]
    }
    pdf(file=plot_file, width=8, height=8, pointsize=10)
        plot(logx, logy, ylab='log(prob > x)', xlab='log(x)', main=plotmain, sub=code)
        abline(coef(xlm))
        legend("topright", lty=c(1, NA),
               legend=c(sprintf("y = %0.2f %0.2f x", coef(xlm)[1], coef(xlm)[2]),
                        sprintf("num snps = %d", nrow(gwas)),
                        sprintf("tail prop = %0.2f", max_prop)))
    dev.off()

    cat(nrow(gwas), "\n", file=donefile)
}
