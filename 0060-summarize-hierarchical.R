#!/usr/bin/env Rscript

library(coda)

source("XX-common.R")

rep <- 1:10

for (drop in c("race", "married", "none")) {

  cat("\n\nOmitted variable(s): ", drop, "\n\n\n")

  prefix <- paste0("out/out-of-sample-hierarchical-0.0000-REP%s/55000-5000-10/5-1-TRUE-1-1-TRUE-TRUE-", drop)

  hosts <- lapply(rep, \(rep) {
    files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
    lapply(files, \(x) readRDS(x)$host)
  })
  print(hosts |> unlist() |> table())

  log_likelihood_samples <- lapply(rep, \(rep) {
    files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
    results <- lapply(files, readRDS)
    log_likelihood_list <- lapply(results, \(x) x$fit$log_likelihood)
    apply(Reduce(cbind, log_likelihood_list), 1, sum)
  })

  doit <- function(log_likelihood) {
    n_eff <- effectiveSize(log_likelihood)
    # Log-Likelihood (from k-fold cross validation)
    me <- qnorm(0.975) * sqrt(var(log_likelihood) / n_eff)
    x <- c(mean(log_likelihood) + c(-1, 0, 1) * me, me)
    names(x) <- c("Lower", "Mean", "Upper", "Margin of Error")
    # print(x)
    x[2]
  }

  totals <- sapply(rep, \(i) doit(log_likelihood_samples[[i]]))
  options(scipen = 9, digits = 16)
  print(mean(totals))
  totals <- transform(totals)
  print(t.test(totals))

}


cat("\n\nOmitted variable(s): ", drop, "\n\n\n")

image_directory <- "out-images"
dir.create(image_directory, showWarnings = FALSE)

files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
results <- lapply(files, readRDS)
k <- length(results)


## Permutation acceptance rate per attempt
p <- mean(sapply(results, \(a) a$fit$permutation_acceptance_rate))
p


## Permutation: proportion of items changed per MCMC iteration
proportion_changed <- mean(sapply(rep, \(rep) {
  files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
  results <- lapply(files, readRDS)
  permutation_list <- lapply(results, \(x) x$fit$permutation)
  mean(sapply(permutation_list, \(permutation) {
    mean(sapply(seq_len(ncol(permutation)-1), \(i) {
      mean(permutation[, i] != permutation[, i + 1])
    }))
  }))
}))
proportion_changed

