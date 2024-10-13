#!/usr/bin/env Rscript

source("XX-common.R")

args <- commandArgs(TRUE)
replication_number <- args[1]

out = "out/ols"
dir.create(out, showWarnings = FALSE)
filename <- file.path(out, paste0(replication_number, ".rds"))

if (file.exists(filename)) {
  cat(paste0(filename, " already exists.\n"))
  q(save = "no")
}

config <- prepare_for_dependent(replication_number, 1, 1, 1, 1, "1994-2020", as.double(0.0))

options(scipen = 9)
is_out_of_sample <- grepl("[0-9]+_[0-9]+", replication_number)

mkData <- function(year, x, returnType = "list") {
  if (nrow(x) == 0) return(NULL)
  attach(x)
  y <- log(AHE)
  X <- cbind(AGE, AGE2 = AGE^2, MALE, WHITE, HISPANIC, UNION, MARRIED, WEEKLYHRS, WEEKLYHRS2 = WEEKLYHRS^2)
  W <- cbind(INTERCEPT = 1, HSGRAD, SOMECOL, BACHELORSPLUS = BACHELORS + ADVDEGREE)   # HSDROPOUT is omitted category
  STATE <- x$KEY
  YEAR <- year
  detach(x)
  if (returnType == "list") {
    list(response = y, global_covariates = X, clustered_covariates = W, item_sizes = table(x$KEY))
  } else {
    as.data.frame(cbind(y, YEAR, STATE, X, W))
  }
}

data <- Reduce(rbind, lapply(seq_along(config$training_data_list), \(i) mkData(i, config$training_data_list[[i]], returnType = "data.frame")))
newdata <- Reduce(rbind, lapply(seq_along(config$validation_data_list), \(i) mkData(i, config$validation_data_list[[i]], returnType = "data.frame")))

fm <- lm(y ~ -1 + as.factor(YEAR):(AGE + AGE2 + MALE + WHITE + HISPANIC + UNION + MARRIED + WEEKLYHRS + WEEKLYHRS2) + as.factor(YEAR):as.factor(STATE):(INTERCEPT + HSGRAD + SOMECOL + BACHELORSPLUS), data = data)
summary(fm)

mean <- predict(fm, newdata, rankdeficient = "NA")
sd <- summary(fm)$sigma

result <- list(missing = sum(is.na(mean)),
               logLike = if (replication_number == "1") NA else sum(dnorm(newdata$y, mean = mean, sd = sd, log = TRUE), na.rm = TRUE),
               time = proc.time())

saveRDS(result, file = filename)
