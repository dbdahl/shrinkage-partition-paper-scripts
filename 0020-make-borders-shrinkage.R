#!/usr/bin/env Rscript

source("XX-common.R")
sessionInfo()

map <- map()

dir.create("shrinkages", showWarnings = FALSE)

map$shrinkage <- 5
map$shrinkage[map$STATENAME %in% c("Maryland", "Delaware", "District of Columbia")] <- 1
map$shrinkage[map$STATENAME %in% c("Montana", "North Dakota", "South Dakota")] <- 1

cat(paste(map$shrinkage, sep = " "), file = "shrinkages/borders.txt")

