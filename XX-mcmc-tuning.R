#!/usr/bin/env Rscript

.before <- ls()

N_ITERATIONS <- 55000
BURNIN       <-  5000
THIN         <-    10

n_replicates <- 10
n_items      <- 51
years        <- 1994:2020

.after <- ls()
.new <- setdiff(.after, c(".before", .before))

args <- commandArgs(TRUE)

if (length(args) == 1 && args[1] == "fish") {
  for (.x in .new) {
    .value <- paste0(eval(parse(text = .x)), collapse = " ")
    cat(paste0("set ", .x," ", .value,"\n"))
  }
}
