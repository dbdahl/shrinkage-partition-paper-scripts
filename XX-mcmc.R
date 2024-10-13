#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
# args <- c("2020", "5500", "500", "10", "fixed", "random", "1", "n", "0", "0.0", "NA", "10_10", "TRUE")
args_original <- args

if (length(args) != 13) {
  cat("usage: YEAR N_ITERATIONS BURNIN THIN PARTITION_DISTRIBUTION PERMUTATION GRIT ANCHOR SHRINKAGE DISCARD_PROPORTION MISSING_ITEM REPLICATION_NUMBER SAVE_SAMPLES")
  q(status = 1)
}

source("XX-common.R")

year <- args_str()
n_iterations <- args_i32()
burnin <- args_i32()
thin <- args_i32()
partition_distribution <- args_str()
permutation <- args_str()
grit <- args_str()
anchor <- args_str()
shrinkage <- args_str()
discard_proportion <- args_str()
missing_item <- args_i32()
replication_number <- args_str()
save_samples <- args_bool()

filename <- file.path(
    if (grepl("[0-9]+_[0-9]+", replication_number)) "out/out-of-sample" else "out/in-sample", year,
    paste(n_iterations, burnin, thin, sep = "-"),
    paste(partition_distribution, grit, anchor, shrinkage, discard_proportion, sep = "-"),
    replication_number,
    paste0(missing_item, ".rds")
)
if (is.na(missing_item)) missing_item <- integer()

if (file.exists(filename)) {
  cat(sprintf("File '%s' already exists.\n", filename))
  q()
}
dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)

config <- config(year, partition_distribution, permutation, grit, anchor, shrinkage, replication_number, as.double(discard_proportion))

mcmc_tuning <- list(
    update_precision_response = TRUE,
    update_global_coefficients = TRUE,
    update_clustering = TRUE,
    update_clustered_coefficients = TRUE,
    n_permutation_updates_per_scan = 10,
    n_items_per_permutation_update = if (config$randomize_permutation) 10 else 0,
    shrinkage_slice_step_size = if (!config$randomize_shrinkage) NULL else {
         a <- config$shrinkage_hyperparameters$shape
         b <- config$shrinkage_hyperparameters$rate
         sqrt(a/b^2)
      },
    grit_slice_step_size = if (!config$randomize_grit) NULL else {
           a <- config$grit_hyperparameters$shape1
           b <- config$grit_hyperparameters$shape2
            sqrt(a*b/((a+b)^2 * (a+b+1)))
      }
    )

mkData <- function(x) {
  if (nrow(x) == 0) return(NULL)
  attach(x)
  y <- log(AHE)
  X <- cbind(AGE, AGE2 = AGE^2, MALE, WHITE, HISPANIC, UNION, MARRIED, WEEKLYHRS, WEEKLYHRS2 = WEEKLYHRS^2)
  W <- cbind(INTERCEPT = 1, HSGRAD, SOMECOL, BACHELORSPLUS = BACHELORS + ADVDEGREE)   # HSDROPOUT is omitted category
  region <- REGION
  detach(x)
  list(response = y, global_covariates = X, clustered_covariates = W, item_sizes = table(x$KEY), region = region)
}

training_data <- mkData(config$training_data)
validation_data <- mkData(config$validation_data)

n_global_coefficients <- ncol(training_data$global_covariates)
n_clustered_coefficients <- ncol(training_data$clustered_covariates)

training_data$region <- NULL
validation_data$region <- NULL
stdev_guess <- 0.361
clustered_coefficients_mean_prior <- matrix(c(1.46, 0.15, 0.24, 0.41), ncol = 1)
clustered_coefficients_precision_prior <- 100 * diag(n_clustered_coefficients)

hyperparameters <- list(precision_response_shape = 1/stdev_guess^2,
                        precision_response_rate = 1,
                        global_coefficients_mean = matrix(rep(0, n_global_coefficients), ncol = 1),
                        global_coefficients_precision = 1 * diag(n_global_coefficients),
                        clustered_coefficients_mean = clustered_coefficients_mean_prior,
                        clustered_coefficients_precision = clustered_coefficients_precision_prior,
                        shrinkage = config$shrinkage_hyperparameters,
                        grit = config$grit_hyperparameters)

global_coefficients_sigma <- solve(hyperparameters$global_coefficients_precision)
clustered_coefficients_sigma <- solve(hyperparameters$clustered_coefficients_precision)

clustering <- local({
    # pd <- config$partition_distribution
    # n_samples <- if (inherits(pd, "CenteredPartition")) 100 else 1
    # x <- samplePartition(pd, n_samples)
    # as.vector(x[n_samples, ])
    config$anchor
})
state <- list(precision_response = rgamma(1, hyperparameters$precision_response_shape, hyperparameters$precision_response_rate),
              global_coefficients = t(mvtnorm::rmvnorm(1, mean = hyperparameters$global_coefficients_mean, sigma = global_coefficients_sigma)),
              clustering = clustering,
              clustered_coefficients = lapply(seq_len(max(clustering)), \(label) {
                  t(mvtnorm::rmvnorm(1, mean = hyperparameters$clustered_coefficients_mean, sigma = clustered_coefficients_sigma))
              }))

compute_log_likelihood <- if (length(missing_item) > 0) {
  if (!is.null(validation_data)) {
    stop("validation_data in the presense of missing items.")
  }
  "missing"
} else if (is.null(validation_data)) {
  "all"
} else {
  "validation"
}

fit <- gourd:::fit(data = training_data, state = state, hyperparameters = hyperparameters,
           partitionDistribution = config$partition_distribution,
           nIterations = n_iterations, burnin = burnin, thin = thin, mcmcTuning = mcmc_tuning,
           missingItems = missing_item, validationData = validation_data,
           save = list(samples = save_samples, logLikelihoodContributions = compute_log_likelihood), progress = FALSE)

out <- list(fit = fit, args = args_original, time = proc.time(), host = Sys.info()[["nodename"]], gourd_version = packageVersion("gourd"))
saveRDS(out, file = filename, compress = FALSE)

