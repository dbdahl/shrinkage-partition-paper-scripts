#!/usr/bin/env Rscript

fm <- readRDS("out/ols/1.rds")$fm
cov_matrix <- vcov(fm)
coefficients <- coef(fm)

nyears <- 27
nstates <- 51

names <- names(coefficients)

result <- vector(nyears, mode = "list")

expected_sqrt_chisq <- function(p) {
  sqrt(2) * gamma((p + 1) / 2) / gamma(p / 2)
}

for (year in seq_len(nyears)) {
  m <- matrix(0, nstates, nstates)
  for (state1 in seq(2, nstates)) {
    for (state2 in seq(1, state1 - 1)) {
      first_set <- which(grepl(sprintf("as[.]factor[(]YEAR[])]%s:as.factor[(]STATE[)]%s:.*", year, state1), names))
      beta1 <- coefficients[first_set]
      include1 <- which(!is.na(beta1))

      second_set <- which(grepl(sprintf("as[.]factor[(]YEAR[])]%s:as.factor[(]STATE[)]%s:.*", year, state2), names))
      beta2 <- coefficients[second_set]
      include2 <- which(!is.na(beta2))

      include <- intersect(include1, include2)
      first_set <- first_set[include]
      second_set <- second_set[include]

      cov_first <- cov_matrix[first_set, first_set]

      # Step 4: Compute the difference vector
      delta <- beta1 - beta2

      # Step 5: Extract the covariance matrices
      cov_second <- cov_matrix[second_set, second_set]
      cov_between <- cov_matrix[first_set, second_set]

      # Step 6: Form the combined covariance matrix
      combined_cov <- cov_first + cov_second - 2 * cov_between

      # Step 7: Solve the system using Cholesky decomposition
      bias <- 0
      while (TRUE) {
        chol_decomp <- tryCatch(chol(combined_cov+bias), error = function(e) NULL)  # Cholesky decomposition
        if (is.null(chol_decomp)) {
          bias <- bias + 0.001
          if (bias > 1) stop("A bridge to far.")
        } else {
          break
        }
      }
      mahalanobis_distance <- if (is.null(chol_decomp)) {browser(); NA} else sqrt(crossprod(forwardsolve(t(chol_decomp), delta))) / expected_sqrt_chisq(length(include))

      # Save the result
      m[state1, state2] <- as.vector(mahalanobis_distance)
    }
  }
  result[[year]] <- m + t(m)
}

saveRDS(result, "out/ols/distances.rds")

