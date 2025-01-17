#!/usr/bin/env Rscript

library(salso)
library(tikzDevice)
library(MASS)
library(gourd)
library(fields)

regions <- scan("anchors/regions.txt")
n_items <- length(regions)
onecluster <- rep(1, n_items)

outfilename <- "0230-summarize-simulation"
outfilenames <- paste0(outfilename, "-", c("vi", "ri"), ".rds")
names(outfilenames) <- c("vi", "ri")

results <- if (all(file.exists(outfilenames))) {
  list(vi = readRDS(outfilenames['vi']), ri = readRDS(outfilenames['ri']))
} else {
  n_seq <- c(100, 400)
  methods <- c(
"sp_jl-2_2-regions-4_1-0.0000",
"cpp_jl_binder-1-regions-198.94-0.0000",
"cpp_jl_binder-1-regions-202.54-0.0000",
"cpp_jl_vi-1-regions-12.43-0.0000",
"cpp_jl_vi-1-regions-18.14-0.0000",
"lsp-1-regions-4.03_1-0.0000",
"lsp-1-regions-7.10_1-0.0000",
"jl-1-regions-0-0.0000",
"sp_jl-2_2-chaos-4_1-0.0000",
"cpp_jl_binder-1-chaos-261.2-0.0000",
"cpp_jl_binder-1-chaos-262.4-0.0000",
"cpp_jl_vi-1-chaos-18.10-0.0000",
"cpp_jl_vi-1-chaos-23.32-0.0000",
"lsp-1-chaos-4.04_1-0.0000",
"lsp-1-chaos-7.10_1-0.0000",
"sp_jl-2_2-chaos-0.0-0.0000",
"sp_jl-2_2-chaos-1.0-0.0000",
"sp_jl-2_2-chaos-2.0-0.0000",
"sp_jl-2_2-chaos-3.0-0.0000",
"sp_jl-2_2-chaos-4.0-0.0000",
"sp_jl-2_2-chaos-5.0-0.0000",
"sp_jl-2_2-chaos-6.0-0.0000",
"sp_jl-2_2-regions-0.0-0.0000",
"sp_jl-2_2-regions-1.0-0.0000",
"sp_jl-2_2-regions-2.0-0.0000",
"sp_jl-2_2-regions-3.0-0.0000",
"sp_jl-2_2-regions-4.0-0.0000",
"sp_jl-2_2-regions-5.0-0.0000",
"sp_jl-2_2-regions-6.0-0.0000"
)
  engine <- function(loss, loss_str) {
    results <- lapply(n_seq, \(n) {
      filenames <- system2("fd", c("-g", sprintf("'d2-sim-%s-*'", n), "out/in-sample/data-clean"), stdout = TRUE)
      results <- methods |> lapply(\(method) {
        ris <- filenames |> sapply(\(x) {
          name <- paste0(x, "/55000-5000-10/", method, "/1/NA.rds")
          if (!file.exists(name)) {
            stop(sprintf("Problem: This file is missing: %s", name))
            return(c(n=0,mean=0,LB=0,UB=0))
          }
          d <- readRDS(paste0(x, "/55000-5000-10/", method, "/1/NA.rds"))
          loss(d$fit$samples$clustering, regions)
        })
        test <- tryCatch(t.test(ris), error = \(x) { list(estimate = mean(ris), conf.int = c(NA, NA)) })
        r <- c(n, test$estimate, test$conf.int)
        names(r) <- c("n", "mean", "LB", "UB")
        r
      })
      names(results) <- methods
      results2 <- Reduce(rbind, results)
      rownames(results2) <- names(results)
      results2
    })
    saveRDS(results, paste0(outfilename, "-", loss_str, ".rds"))
    results
  }
  list(vi = engine(VI, "vi"), ri = engine(RI, "ri"))
}

results

## Shrinkage and grit
image_directory <- "out-images"
dir.create(image_directory, showWarnings = FALSE)

x <- readRDS(paste0("out/in-sample/data-clean/d2-sim-400-1.rds/55000-5000-10/sp_jl-2_2-regions-4_1-0.0000/1/NA.rds"))
y <- readRDS(paste0("out/in-sample/data-clean/d2-sim-400-1.rds/55000-5000-10/sp_jl-2_2-chaos-4_1-0.0000/1/NA.rds"))

all_shrinkages_x <- x$fit$samples$shrinkage[, 1]
all_shrinkages_y <- y$fit$samples$shrinkage[, 1]
shrinkage_alpha_x <- 4
shrinkage_alpha_y <- 4
shrinkage_beta <- 1
shrinkage_mean_prior_x <- shrinkage_alpha_x / shrinkage_beta
shrinkage_mean_prior_y <- shrinkage_alpha_y / shrinkage_beta
shrinkage_mean_posterior_x <- mean(all_shrinkages_x)
shrinkage_mean_posterior_y <- mean(all_shrinkages_y)
cat("Prior mean: ", shrinkage_mean_prior_x, "\n", sep = "")
cat("Posterior mean: ", shrinkage_mean_posterior_x, "\n", sep = "")
cat("Prior mean: ", shrinkage_mean_prior_y, "\n", sep = "")
cat("Posterior mean: ", shrinkage_mean_posterior_y, "\n", sep = "")

all_grits_x <- x$fit$samples$grit
all_grits_y <- y$fit$samples$grit
grit_alpha <- 2
grit_beta <- 2
grit_mean_prior <- grit_alpha / (grit_alpha + grit_beta)
grit_mean_posterior_x <- mean(all_grits_x)
grit_mean_posterior_y <- mean(all_grits_y)
cat("Prior mean: ", grit_mean_prior, "\n", sep = "")
cat("Posterior mean: ", grit_mean_posterior_x, "\n", sep = "")
cat("Posterior mean: ", grit_mean_posterior_y, "\n", sep = "")
  
kde_x <- kde2d(all_shrinkages_x, all_grits_x, n = 100)
kde_y <- kde2d(all_shrinkages_y, all_grits_y, n = 100)

set.seed(34534512)
n_mc_samples <- 1000
shrinkage_n <- 100
grit_n <- 100
#   n_mc_samples <- 100
#   shrinkage_n <- 25
#   grit_n <- 25

out <- summarize_prior_on_shrinkage_and_grit(regions, n_mc_samples = n_mc_samples, shrinkage_n = shrinkage_n, shrinkage_shape = shrinkage_alpha_x, shrinkage_rate = shrinkage_beta, grit_n = grit_n, grit_shape1 = grit_alpha, grit_shape2 = grit_beta, use_crp = FALSE, concentration = 1.0, domain_specification = list(shrinkage_lim = c(0, 10), grit_lim = c(0.0, 1.0)))

shrinkage_grit_plot <- function(background) {
  png(file.path(image_directory, sprintf("shrinkage_grit_%s.png", background)), width = 1500, height = 1000, pointsize = 64) # pointsize = 48
  if (background == "posterior") {
    par(mar = c(4, 4, 0.2, 5) + 0.1, las = 1, family = "serif", xaxs = "i", yaxs = "i", plt = c(0.170, 0.82, 0.27, 0.97696))
    print(par(no.readonly = TRUE)$plt)
    plot(NA, xlab = "", ylab = "", type = "n",
      xlim = range(out$shrinkage) + c(-1,1)*(diff(range(out$shrinkage))*0.005),
      ylim = range(out$grit) + c(-1,1)*(diff(range(out$grit))*0.005))
    mtext("Shrinkage", 1, line = 2.5, las = 0)
    mtext("Grit", 2, line = 2.5, las = 0)
    contour(kde_x, add = TRUE, lwd = 6, lty = 1, drawlabels = FALSE, col = "black")
    contour(kde_y, add = TRUE, lwd = 6, lty = 5, drawlabels = FALSE, col = "black")
    arrows(
      shrinkage_mean_prior_x, grit_mean_prior,
      shrinkage_mean_posterior_x, grit_mean_posterior_x,
      lwd = 14, cex = 1.6, pch = 21, col = "black", length = 0.4
      )
    arrows(
      shrinkage_mean_prior_y, grit_mean_prior,
      shrinkage_mean_posterior_y, grit_mean_posterior_y,
      lwd = 14, cex = 1.6, pch = 21, col = "black", length = 0.4
      )
  } else {
    par(mar = c(4, 4, 0.2, 2) + 0.1, las = 1, family = "serif", xaxs = "i", yaxs = "i")
    if (background == "expected_rand_index") {
      by <- 0.008
      breaks <- c(0, seq(0.50 - by / 2, 1.0 + by / 2, by = by))
      zlim <- c(0, 1)
      image.plot(out$shrinkage, out$grit, out[[background]], xlab = "", ylab = "", legend.width = 0.7, legend.mar = 4, bigplot = c(0.170, 0.82, 0.27, 0.97696), smallplot = c(0.84, 0.869, 0.295, 0.943856), zlim = zlim, breaks = breaks, nlevel = length(breaks)-1)
    } else {
      zlim <- c(0, 2.7)
      image.plot(out$shrinkage, out$grit, out[[background]], xlab = "", ylab = "", legend.width = 0.7, legend.mar = 4, bigplot = c(0.170, 0.82, 0.27, 0.97696), smallplot = c(0.84, 0.869, 0.295, 0.943856), zlim = zlim)
    }
    mtext("Shrinkage", 1, line = 2.5, las = 0)
    mtext("Grit", 2, line = 2.5, las = 0)
    contour(out$shrinkage, out$grit, exp(out$log_density), add = TRUE, lwd = 6, lty = 1, drawlabels = FALSE, col = "white")
  }
  dev.off()
}

shrinkage_grit_plot("expected_rand_index")
shrinkage_grit_plot("expected_entropy")
shrinkage_grit_plot("posterior")

tikz(file.path(image_directory, "simulation_1.tex"), width = 4, height = 2.7)
shrink_seq <- 0:6
chaos <- sprintf("sp_jl-2_2-1-%s.0-0.0000", shrink_seq)
regions <- sprintf("sp_jl-2_2-regions-%s.0-0.0000", shrink_seq)
par(mar = c(4, 4, 0.4, 0.2))
plot(NA, xlim = c(0, max(shrink_seq)), ylim = c(0.6, 1.0), xlab = "Shrinkage $\\omega$", ylab = "E(Rand Index $\\mid$ Data)",
     axes = FALSE, xaxs = "i", yaxs = "i")
axis(2, las = 1)
axis(1, las = 1)
box()
tmp <- results$ri[[1]][row.names(results$ri[[1]]) %in% sprintf("sp_jl-2_2-regions-%0.1f-0.0000", shrink_seq), "mean"]
points(shrink_seq, tmp, cex = 1, lty = 3)
lines(shrink_seq, tmp, cex = 1, lty = 3)
tmp <- results$ri[[2]][row.names(results$ri[[2]]) %in% sprintf("sp_jl-2_2-regions-%0.1f-0.0000", shrink_seq), "mean"]
points(shrink_seq, tmp, cex = 1, lty = 1)
lines(shrink_seq, tmp, cex = 1, lty = 1)
tmp <- results$ri[[1]][row.names(results$ri[[1]]) %in% sprintf("sp_jl-2_2-chaos-%0.1f-0.0000", shrink_seq), "mean"]
points(shrink_seq, tmp, cex = 1, lty = 3, pch = 8)
lines(shrink_seq, tmp, cex = 1, lty = 3)
tmp <- results$ri[[2]][row.names(results$ri[[2]]) %in% sprintf("sp_jl-2_2-chaos-%0.1f-0.0000", shrink_seq), "mean"]
points(shrink_seq, tmp, cex = 1, lty = 1, pch = 8)
lines(shrink_seq, tmp, cex = 1, lty = 1)
legend(0.10, 1.00, c("$m=100$, $\\mu$ is regions", "$m=100$, $\\mu$ is shuffled", "$m=400$, $\\mu$ is regions", "$m=400$, $\\mu$ is shuffled"), pch=c(1,8,1,8), lty=c(3,3,1,1), cex=0.85, bty="n")
dev.off()

