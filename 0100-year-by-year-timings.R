#!/usr/bin/env Rscript

template <- "out/in-sample/%s/55000-5000-10/%s/1/NA.rds"
methods <-  c("cpp_crp_binder-1-regions-500-0.0000", "cpp_crp_vi-1-regions-30-0.0000", "crp-1-regions-0-0.0000", "fixed-1-regions-0-0.0000", "lsp-1-regions-3990.81-0.0000", "sp_crp-0.02-regions-5-0.0000", "sp_crp-0.02-regions-borders-0.0000")
years <- 1994:2020

results <- lapply(methods, \(m) {
  z <- lapply(years, \(y) {
    x <- readRDS(sprintf(template, y, m))
    list(host = x$host, time = sum(x$time[c(1, 4)]), n_clusters = mean(apply(x$fit$samples$clustering, 1, \(x) x |> unique() |> length())))
  })
  list(host = table(sapply(z, \(zz) zz$host)), time = sum(sapply(z, \(zz) zz$time)), n_clusters = mean(sapply(z, \(zz) zz$n_clusters)))
})

sapply(results, \(x) x$host)

times <- sapply(results, \(x) x$time)
names(times) <- methods
times / 60  # Total time in minutes

n_clusters <- sapply(results, \(x) x$n_clusters)
names(n_clusters) <- methods
n_clusters

