#!/usr/bin/env Rscript

options(scipen = 9, digits = 16)

files <- system2("fd", c("-I", "-g", "'*_*.rds'", "out/ols"), stdout = TRUE)
x <- lapply(files, readRDS)

n_observations_excluded <- x |> sapply(\(y) y$missing) |> sum()
ols_fit <- x |> sapply(\(y) y$logLike) |> sum()
ols_fit

source("XX-common.R")
transform(ols_fit, 139555 - n_observations_excluded)

