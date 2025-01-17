#!/usr/bin/env Rscript

library(salso)
library(fields)
library(dynamicTreeCut)

distances <- readRDS("out/ols/distances.rds")

trees <- lapply(distances, \(dist) hclust(as.dist(dist), method = "complete"))

clusterings <- vector(mode = "list", length = length(trees))
for (i in seq_along(trees)) {
  clusterings[[i]] <- cutreeDynamic(dendro = trees[[i]], distM = distances[[i]])
}

m <- matrix(0.0, nrow = length(clusterings), ncol = length(clusterings))
for (i in seq_along(clusterings)) {
  for (j in seq_along(clusterings)) {
    m[i,j] <- RI(clusterings[[i]], clusterings[[j]])
  }
}

image_directory <- "out-images"
filename <- file.path(image_directory, "ri_matrix_hclust.png")
png(filename, width = 1200, height = 1200, pointsize = 48)
par(mar = c(2.2, 3.1, 0.1, 0.7), family = "serif")
image(m[, rev(seq_len(nrow(m)))], axes = FALSE)
box()
axis(1, at = c(1, 6, 11, 16, 21, 26) / 26, label = c(1995, 2000, 2005, 2010, 2015, 2020), las = 1)
axis(2, at = c(0, 5, 10, 15, 20, 25) / 26, label = rev(c(1995, 2000, 2005, 2010, 2015, 2020)), las = 1)
dev.off()
