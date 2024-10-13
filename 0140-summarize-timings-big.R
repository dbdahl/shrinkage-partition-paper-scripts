#!/usr/bin/env Rscript

library(txtplot)

extract <- function(files) {
  data <- lapply(files, \(y) readRDS(y))
  host <- sapply(data, \(y) y$host)
  time <- sapply(data, \(y) y$time[['elapsed']])
  size <- as.integer(sub(".*-mega-([0-9]+).*", "\\1", files))
  n_clusters <- sapply(data, \(data) mean(apply(data$fit$samples$clustering, 1, \(x) length(unique(x)))))
  list(host = host, time = time, size = size, n_clusters)
}

files <- system2("fd",c("-I","-p","'out/in-sample/mega-.*/55000-5000-10/crp-1-regions-mega-.*-0-0.0000/.*/NA.rds'"), stdout = TRUE)
crp <- extract(files)

files <- system2("fd",c("-I","-p","'out/in-sample/mega-.*/55000-5000-10/sp_crp-0.02-regions-mega-.*-5-0.0000/.*/NA.rds'"), stdout = TRUE)
sp <- extract(files)

table(crp$host)
table(sp$host)

with(crp, txtplot(size, time))
with(sp, txtplot(size, time))

both <- merge(as.data.frame(crp), as.data.frame(sp), by = "size")[, c(1,3,4,6,7)]
names(both) <- c("size","crp_time","crp_n_clusters","sp_time","sp_n_clusters")
both$ratio <- both$sp_time / both$crp_time
both$crp_time_in_hours <- both$crp_time / 60^2
both$sp_time_in_hours <- both$sp_time / 60^2
both

both$n <- 51*both$size

plot(both$n, both$sp_time, type = "l")
lines(both$n, both$crp_time)

plot(both$n, log(both$sp_time), type = "l")
lines(both$n, log(both$crp_time))

plot(both$n, sqrt(both$sp_time), type = "l")
lines(both$n, sqrt(both$crp_time))

fm <- lm(sqrt(crp_time) ~ n, data = both)
summary(fm)
hat <- predict(fm)
lines(both$n, hat)

fm <- lm(sqrt(sp_time) ~ n, data = both)
summary(fm)
hat <- predict(fm)
lines(both$n, hat)
