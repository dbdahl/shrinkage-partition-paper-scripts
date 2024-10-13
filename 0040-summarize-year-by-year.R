#!/usr/bin/env Rscript

source("XX-common.R")

if (!file.exists("0040-summarize-year-by-year.RData")) {

  data <- read_data()
  n_observations <- nrow(data)

  filename <- system2("fd", c("\\.rds", "out/out-of-sample"), stdout = TRUE)
  distribution <- basename(dirname(dirname(filename)))
  fold_info <- matrix(as.integer(unlist(strsplit(basename(dirname(filename)), "_"))), ncol = 2, byrow = TRUE)
  fold <- fold_info[, 1]
  n_folds <- unique(as.integer(fold_info[, 2]))
  if (length(n_folds) != 1) stop("Inconsistent number of folds.")
  year <- as.integer(basename(dirname(dirname(dirname(dirname(filename))))))
  unique_distribution <- unique(distribution)
  e <- lapply(unique_distribution, \(distr) {
    m <- Reduce(`+`, lapply(filename[distribution == distr], \(file) {
      readRDS(file)$fit$logLikelihoodContribution
    }))
    discard_proportion <- tools::file_ext(distr)
    list(apply(m, 1, sum), apply(m, 2, mean), discard_proportion = as.integer(discard_proportion)/10^(nchar(discard_proportion)))
  })

  discard_proportion <- sapply(e, \(ee) ee$discard_proportion)
  ci <- t(sapply(lapply(e, \(ee) ee[[1]]), \(x) {
    n <- coda::effectiveSize(x)
    result <- mean(x) + c( -1, 0, 1) * qnorm(0.975) * sqrt(var(x) / n)
    names(result) <- c("lower", "mean", "upper")
    result
  }))
  f <- data.frame(distribution = unique_distribution, ci, n_observations = (1 - discard_proportion) * n_observations)
  f$margin_of_error <- f$upper - f$mean
  ft <- transform(f[, c(2,3,4)])
  names(ft) <- paste0(names(ft), "_")
  f <- cbind(f, ft)
  f$margin_of_error_ <- f$upper_ - f$mean_
  f <- f[order(f$mean_, decreasing = TRUE), ]
  g <- f

  save.image("0040-summarize-year-by-year.RData")

} else {

  load("0040-summarize-year-by-year.RData")
  
}

options(width = 300, scipen = 9, digits = 16)
g
max(g$margin_of_error_)

h <- g[, "mean_", drop = FALSE]
row.names(h) <- g[, 1]
round(h)

featured <- c("crp-1-regions-0-0.0000", "fixed-1-regions-0-0.0000", "cpp_crp_vi-1-regions-30-0.0000", "sp_crp-0.02-regions-5-0.0000", "sp_crp-0.02-regions-borders-0.0000", "lsp-1-regions-3990.81-0.0000", "cpp_crp_binder-1-regions-500-0.0000")
round(h[row.names(h) %in% featured, , drop = FALSE])


n_items <- 51
w <- cbind(seq_len(n_items), Reduce(cbind, lapply(e, \(ee) ee[[2]])))
colnames(w) <- c("KEY", unique_distribution)

w0 <- w

colnames(w)
w <- w[, c("KEY", featured)]

options(scipen = 2)
map <- map()
map$STATECENSUS <- NULL
m <- merge(w, map)
m

m2 <- aggregate(m[,setdiff(names(m),c("KEY","STATENAME","REGION","DIVISION"))], list(REGION = m[,"REGION"]), FUN=sum)
m2

apply(m2, 2, sum)


