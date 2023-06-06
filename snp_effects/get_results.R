#!/usr/bin/env Rscript

## Original run with tail prob 25%
results <- data.frame(code=setdiff(list.dirs(".", recursive=FALSE, full.names=FALSE), ".git"), stringsAsFactors=FALSE)
coefs <- t(sapply(results$code, function (x) { f <- file.path(x, "xlm.RData"); if (file.exists(f)) { load(f); coef(xlm) } else { c(NA, NA) } } ))
results$intercept <- coefs[,1]
results$slope <- coefs[,2]
results$nsnps <- sapply(results$code, function (x) { f <- file.path(x, "done"); if (file.exists(f)) { scan(f) } else { NA } } )
write.table(results, file="snp_results_25.tsv")

## Next run with tail prob 10%
results <- data.frame(code=setdiff(list.dirs(".", recursive=FALSE, full.names=FALSE), ".git"), stringsAsFactors=FALSE)

max_prop <- 0.1
filebase <- sprintf("tail_%d_", floor(100*max_prop))
basename <- file.path(results$code, filebase)
donefile <- paste0(basename, "done")
xlm_file <- paste0(basename, "xlm.RData")

coefs <- t(sapply(xlm_file, function (f) { if (file.exists(f)) { load(f); coef(xlm) } else { c(NA, NA) } } ))
results$intercept <- coefs[,1]
results$slope <- coefs[,2]
results$nsnps <- sapply(donefile, function (f) { if (file.exists(f)) { scan(f) } else { NA } } )
write.table(results, file="snp_results_10.tsv")


