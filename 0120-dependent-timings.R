#!/usr/bin/env Rscript

library(fields)
library(MASS)
library(gourd)
library(salso)  # Requires version 0.3.38 or higher.
# library(tikzDevice)

for (model in c("hierarchical", "temporal")) {
  for (years in c("1994-2002", "1994-2011", "1994-2020", "1994-2020-2011", "1994-2020-2002", "1994-2020-1994")) {
    files <- system2("fd",c("-I","-g","'*.rds'",paste0("out/in-sample-", model, "-", years)), stdout = TRUE)
    files <- files[grepl(".*/[1-5].rds", files)] ## Only take the first 5 for consistency
    files |> sapply(\(x) readRDS(x)$time[['user.self']]) |> mean() |> print()
    files |> sapply(\(x) readRDS(x)$host) |> table() |> print()
    time <- files |> sapply(\(x) readRDS(x)$time[['elapsed']]) |> mean() # Use this one... but should be similar to above.
    cat(paste0("Years: ", years, "; Model: ", model, "; Mean time: ", time, "\n"))
  }
}

# How time is apportioned
model <- "hierarchical"
years <- "1994-2020"
files <- system2("fd",c("-I","-g","'*.rds'",paste0("out/in-sample-", model, "-", years)), stdout = TRUE)
files <- files[grepl(".*/[1-5].rds", files)] ## Only take the first 5 for consistency
y <- apply(Reduce(rbind, lapply(files, \(x) readRDS(x)$fit$wall_time)), 2, mean)
y / sum(y)

# How time is apportioned
model <- "temporal"
years <- "1994-2020"
files <- system2("fd",c("-I","-g","'*.rds'",paste0("out/in-sample-", model, "-", years)), stdout = TRUE)
files <- files[grepl(".*/[1-5].rds", files)] ## Only take the first 5 for consistency
y <- apply(Reduce(rbind, lapply(files, \(x) readRDS(x)$fit$wall_time)), 2, mean)
y / sum(y)

year <- "2020"
image_directory <- "out-images"
dir.create(image_directory, showWarnings = FALSE)


## Hierarchical visualizations
model <- "hierarchical"
files <- system2("fd",c("-I","-g","'*.rds'",paste0("out/in-sample-", model, "-1994-", year)), stdout = TRUE)
x <- Reduce(rbind, lapply(files, \(x) t(readRDS(x)$fit$anchor)))
est <- salso(x)
est
table(est)
summ <- summary(est)
# tikz(file.path(image_directory, "anchor.tex"), width = 6, height = 6)
png(file.path(image_directory, "anchor.png"), width = 1200, height = 1200, pointsize = 48)
par(family = "serif")
plot(summ, showLabels = TRUE, cexAdjustment = 0.5)
dev.off()


## Temporal visualizations
model <- "temporal"
source("./XX-common.R")
map <- map()
source("./XX-mcmc-tuning.R", echo = TRUE)

files <- system2("fd",c("-I","-g","'*.rds'",paste0("out/in-sample-", model, "-1994-", year)), stdout = TRUE)
samples <- lapply(files, \(x) readRDS(x)$fit$unit_partitions |> lapply(\(y) t(y)))
samples <- lapply(seq_along(years), \(t) Reduce(rbind, lapply(samples, \(x) x[[t]])))

samples_crp <- lapply(years, \(year) {
  dir <- sprintf("out/in-sample/%s/%s-%s-%s/crp-1-regions-0-0.0000", year, N_ITERATIONS, BURNIN, THIN)
  readRDS(sprintf("%s/1/NA.rds", dir))$fit$samples$clustering
  #dir <- sprintf("out/%s/%s-%s-%s/crp-1-regions-0-0.0000", year, N_ITERATIONS, BURNIN, THIN)
  #Reduce(rbind, lapply(lapply(seq_len(n_replicates), \(i) readRDS(sprintf("%s/%s_%s/NA.rds", dir, i, n_replicates))), \(x) x$fit$samples$clustering))
})

nRuns <- 16
estimates_all <- lapply(samples, \(x) salso(x, nRuns = nRuns))
estimates <- Reduce(rbind, estimates_all)

estimates_crp_all <- lapply(samples_crp, \(x) salso(x, nRuns = nRuns))
estimates_crp <- Reduce(rbind, estimates_crp_all)

colors <- hcl.colors(5, "YlOrRd", rev = TRUE)[c(1,5,3,2,4)]
partition_over_time_plot <-  function(x) {
  # par(mar = c(2.2, 0.1, 0.1, 0.63))  
  par(mar = c(2.2, 0.2, 0.1, 0.2))
  swap <- function(x, i, j) {
    x[x==i] <- -1
    x[x==j] <- i
    x[x==-1] <- j
    x
  }
  x[2,] <- swap(x[2,], 4, 5)
  x[7,] <- swap(x[7,], 3, 4)
  x[8,] <- swap(x[8,], 3, 4)
  x[8,] <- swap(x[8,], 4, 5)
  x[10,] <- swap(x[10,], 3, 4)
  x[10,] <- swap(x[10,], 4, 5)
  x[18,] <- swap(x[18,], 3, 4)
  image(x, col = colors, axes = FALSE)
  box()
  grid <- \(len) seq(0, 1, length = len) + 0.5 / len
  abline(v = grid(27))
  abline(h = grid(51))
  axis(1, at = c(1, 6, 11, 16, 21) / 26, label = c(1995, 2000, 2005, 2010, 2015), las = 1)
}

# tikz(file.path(image_directory, "partitions_over_time.tex"), width = 8, height = 3)
png(file.path(image_directory, "partitions_over_time.png"), width = 1600, height = 600,pointsize = 48)
layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))
partition_over_time_plot(estimates_crp)
partition_over_time_plot(estimates)
dev.off()


pseudo_correlation_plot <- function(filename, samples, loss = RI, max_minus = FALSE, n = 1000) {
  n_years <- length(years)
  lm <- matrix(0.0, nrow = n_years, ncol = n_years)
  for (i in seq_len(n_years)) {
    for (j in seq(i, n_years)) {
      ni <- nrow(samples[[i]])
      nj <- nrow(samples[[j]])
      lm[i, j] <- mean(loss(samples[[i]][sample(ni, n), ], samples[[j]][sample(nj, n), ]))
    }
  }
  lm2 <- t(lm) + lm - diag(diag(lm))   # Symmetrize
  lm4 <- if (max_minus) {
    lm3 <- max(lm2) - lm2 # Loss to pseudo covariance
    cov2cor(lm3)                  # Pseudo covariance to pseudo correlation
  } else {
    lm2
  }
  # tikz(filename, width = 6, height = 6)
  png(filename, width = 1200, height = 1200, pointsize = 48)
  par(mar = c(2.2, 3.1, 0.1, 0.7), family = "serif")
  image(lm4[, rev(seq_len(nrow(lm4)))], axes = FALSE)
  box()
  axis(1, at = c(1, 6, 11, 16, 21, 26) / 26, label = c(1995, 2000, 2005, 2010, 2015, 2020), las = 1)
  axis(2, at = c(0, 5, 10, 15, 20, 25) / 26, label = rev(c(1995, 2000, 2005, 2010, 2015, 2020)), las = 1)
  dev.off()
}

n_samples <- 500
pseudo_correlation_plot(file.path(image_directory, "ri_matrix_sp.png"), samples, n = n_samples, loss = RI)
pseudo_correlation_plot(file.path(image_directory, "ri_matrix_crp.png"), samples_crp, n = n_samples, loss = RI)

