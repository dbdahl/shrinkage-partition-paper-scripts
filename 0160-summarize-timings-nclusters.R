#!/usr/bin/env Rscript

concentration <- c("0.01", "0.1", "1", "5", "10")
sp_crp <- sprintf("sp_crp%s-0.5-regions-5-0.0000", concentration)

files <- function(x) {
  lapply(x, \(d) system2("fd",c("-I","-p", sprintf("'out/in-sample/1994/55000-5000-10/%s/.*/NA.rds'", d)), stdout = TRUE))
}

getdata <- function(files) {
  as.data.frame(Reduce(rbind, lapply(files, \(x) t(sapply(x, \(x) {
      z <- readRDS(x)
      mean <- z$fit$samples$clustering |> apply(1, \(y) length(unique(y))) |> mean()
      time <- z$time['elapsed']
      c(n_clusters = mean, time)
  })))))
}

sp_crp_results <- getdata(files(sp_crp))
sp_crp_results

with(sp_crp_results, {
    plot(n_clusters, elapsed)
    fm <- lm(elapsed ~ n_clusters)
    print(summary(fm))
    hat <- predict(fm)
    lines(n_clusters, hat)
})



