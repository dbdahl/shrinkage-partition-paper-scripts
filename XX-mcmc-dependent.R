#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
#  args <- c("hierarchical", "1994-2020", "5", "1", "TRUE", "2", "2", "TRUE", "TRUE", "race", "550", "50", "1", "10_10", "10")
args_original <- args

if (length(args) != 16) {
  cat("usage: [ hierarchical | temporal ] YEARS DISCARD_PROPORTION SHRINKAGE_SHAPE SHRINKAGE_RATE UPDATE_SHRINKAGE GRIT_SHAPE1 GRIT_SHAPE2 UPDATE_GRIT IDIOSYNCRATIC_PERMUTATION DROP N_ITERATIONS BURNIN THINNING REPLICATION_NUMBER RUN")
  q(status = 1)
}

source("XX-common.R")

model_name <- args_str()
years <- args_str()
discard_proportion <- args_str()
shrinkage_shape <- args_f64()
shrinkage_rate <- args_f64()
update_shrinkage <- args_bool()
grit_shape1 <- args_f64()
grit_shape2 <- args_f64()
update_grit <- args_bool()
idiosyncratic_permutation <- args_bool()
drop <- args_str()
n_iterations <- args_i32()
burnin <- args_i32()
thinning <- args_i32()
replication_number <- args_str()
run_number <- args_i32()

anchor_concentration = 1.0    # Fixed
baseline_concentration = 1.0  # Fixed

config <- prepare_for_dependent(replication_number, shrinkage_shape, shrinkage_rate, grit_shape1, grit_shape2, years, as.double(discard_proportion))

sha <- config$shrinkage_hyperparameters$shape
shb <- config$shrinkage_hyperparameters$rate
sga <- config$grit_hyperparameters$shape1
sgb <- config$grit_hyperparameters$shape2

mcmc_tuning <- list(
    update_precision_response = TRUE,
    update_global_coefficients = TRUE,
    update_clustering = TRUE,
    update_clustered_coefficients = TRUE,
    n_permutation_updates_per_scan = if ( idiosyncratic_permutation ) 10 else 200,
    n_items_per_permutation_update = 5,
    shrinkage_slice_step_size = if (!update_shrinkage) NULL else sqrt(sha/shb^2),
    grit_slice_step_size = if (!update_grit) NULL else sqrt(sga*sgb/((sga+sgb)^2 * (sga+sgb+1)))
    )

options(scipen = 9)
is_out_of_sample <- grepl("[0-9]+_[0-9]+", replication_number)
filename <- file.path(
    paste0(if (is_out_of_sample) "out/out-of-sample-" else "out/in-sample-", model_name, if (is_out_of_sample) paste0("-", discard_proportion, "-REP", run_number) else paste0("-", years)), paste0(n_iterations, "-", burnin, "-", thinning), paste0(sha, "-", shb, "-", update_shrinkage, "-", sga, "-", sgb, "-", update_grit, "-", idiosyncratic_permutation, "-", drop), paste0(replication_number, ".rds")
)

if (file.exists(filename)) {
  cat(sprintf("File '%s' already exists.\n", filename))
  q()
}
dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)

mkData <- function(x) {
  if (nrow(x) == 0) return(NULL)
  attach(x)
  y <- log(AHE)
  X <- if (drop == "married") {
    cbind(AGE, AGE2 = AGE^2, MALE, WHITE, HISPANIC, UNION,          WEEKLYHRS, WEEKLYHRS2 = WEEKLYHRS^2)
  } else if (drop == "race") {
    cbind(AGE, AGE2 = AGE^2, MALE,                  UNION, MARRIED, WEEKLYHRS, WEEKLYHRS2 = WEEKLYHRS^2)
  } else if (drop == "none") {
    cbind(AGE, AGE2 = AGE^2, MALE, WHITE, HISPANIC, UNION, MARRIED, WEEKLYHRS, WEEKLYHRS2 = WEEKLYHRS^2)
  } else {
    stop("Unrecognized 'DROP' argument.")
  }
  W <- cbind(INTERCEPT = 1, HSGRAD, SOMECOL, BACHELORSPLUS = BACHELORS + ADVDEGREE)   # HSDROPOUT is omitted category
  detach(x)
  list(response = y, global_covariates = X, clustered_covariates = W, item_sizes = table(x$KEY))
}

training_data_list <- lapply(config$training_data_list, \(training_data) mkData(training_data))
validation_data_list <- if (is.null(config$validation_data_list)) {
  NULL
} else {
  lapply(config$validation_data_list, \(validation_data) mkData(validation_data))
}

n_global_coefficients <- ncol(training_data_list[[1]]$global_covariates)
n_clustered_coefficients <- ncol(training_data_list[[1]]$clustered_covariates)

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

clustering <- sample(1:5, length(training_data_list[[1]]$item_sizes), replace = TRUE)

state <- list(precision_response = rgamma(1, hyperparameters$precision_response_shape, hyperparameters$precision_response_rate),
              global_coefficients = t(mvtnorm::rmvnorm(1, mean = hyperparameters$global_coefficients_mean, sigma = global_coefficients_sigma)),
              clustering = clustering,
              clustered_coefficients = lapply(seq_len(max(clustering)), \(label) {
                  t(mvtnorm::rmvnorm(1, mean = hyperparameters$clustered_coefficients_mean, sigma = clustered_coefficients_sigma))
              }))

unit_list <- lapply(training_data_list, \(training_data) {
  list(training_data = training_data, state = state, hyperparameters = hyperparameters)
})

global_mcmc_tuning <- list(n_iterations = n_iterations, burnin = burnin, thinning = thinning, update_anchor = TRUE, idiosyncratic_permutation = idiosyncratic_permutation)

fit <- gourd:::fit_dependent(model_name, unit_list, config$anchor_anchor, config$anchor_shrinkage, config$anchor_shrinkage_reference, anchor_concentration, baseline_concentration, hyperparameters, mcmc_tuning, global_mcmc_tuning, validation_data_list)

out <- list(fit = fit, args = args_original, time = proc.time(), host = Sys.info()[["nodename"]], gourd_version = packageVersion("gourd"))
saveRDS(out, file = filename, compress = FALSE)

