#!/usr/bin/env Rscript

library(salso)  # Requires version 0.3.38 or higher.
# library(tikzDevice)
library(coda)
library(MASS)
library(fields)

source("XX-common.R")

rep <- 1:10

log_likelihood_samples <- function(rep, prefix) {
  lapply(rep, \(rep) {
    files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
    results <- lapply(files, readRDS)
    log_likelihood_list <- lapply(results, \(x) x$fit$log_likelihood)
    apply(Reduce(cbind, log_likelihood_list), 1, sum)
  })
}

doit <- function(log_likelihood) {
  n_eff <- effectiveSize(log_likelihood)
  # Log-Likelihood (from k-fold cross validation)
  me <- qnorm(0.975) * sqrt(var(log_likelihood) / n_eff)
  x <- c(mean(log_likelihood) + c(-1, 0, 1) * me, me)
  names(x) <- c("Lower", "Mean", "Upper", "Margin of Error")
  # print(x)
  x[2]
}


for (drop in c("race", "married", "none")) {

  cat("\n\nOmitted variable(s): ", drop, "\n\n\n")

  prefix <- paste0("out/out-of-sample-temporal-0.0000-REP%s/55000-5000-10/5-1-TRUE-1-1-TRUE-TRUE-", drop)

  hosts <- lapply(rep, \(rep) {
    files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
    lapply(files, \(x) readRDS(x)$host)
  })
  print(hosts |> unlist() |> table())

  totals <- sapply(rep, \(i) doit(log_likelihood_samples(rep, prefix)[[i]]))
  options(scipen = 9, digits = 16)
  print(mean(totals))
  totals <- transform(totals)
  print(t.test(totals))

}

files <- system2("fd", c("-I", "-g", "'*.rds'", sprintf(prefix, rep)), stdout = TRUE)
all <- lapply(files, readRDS)

## Permutation acceptance rate per attempt
p <- mean(sapply(all, \(a) a$fit$permutation_acceptance_rate))
p


drop <- "none"

for (discard_proportion in c("0.1000", "0.2000", "0.3000", "0.4000", "0.5000")) {

  cat("\n\nDiscard proportion: ", discard_proportion, "\n\n\n")

  prefix <- paste0("out/out-of-sample-temporal-", discard_proportion, "-REP%s/55000-5000-10/5-1-TRUE-1-1-TRUE-TRUE-", drop)

  totals <- sapply(rep, \(i) doit(log_likelihood_samples(rep, prefix)[[i]]))
  options(scipen = 9, digits = 16)
  print(mean(totals))
  totals <- transform(totals)
  print(t.test(totals))

}

