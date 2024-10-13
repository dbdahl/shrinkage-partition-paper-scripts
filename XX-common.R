library(gourd)

options(warn = 2)

args_str <- function() {
  x <- args[1]
  args <<- args[-1]
  x
}

args_i32 <- function() {
  x <- args_str()
  if (x == "NA") NA_integer_ else as.integer(x)
}

args_f64 <- function() {
  x <- args_str()
  if (x == "NA") NA_real_ else as.double(x)
}

args_bool <- function() {
  x <- args_str()
  if (x == "NA") NA else as.logical(x)
}

process_str <- function(str) {
  if (length(str) != 1) stop("Unexpected length for argument.")
  if (endsWith(str, ".R")) {
    list(TRUE, source(str)$value)
  } else if (endsWith(str, ".rds")) {
    list(TRUE, readRDS(str))
  } else if (file.exists(str)) {
    list(TRUE, scan(str, what = numeric()))
  } else {
    FALSE
  }
}

read_data <- function(year = NA) {
  all <- if (!is.na(year) && isTRUE(endsWith(year, ".rds"))) {
    readRDS(year)
  } else {
    all <- readRDS("data-clean/d2.rds")
    value <- sort(unique(all$STATECENSUS))
    map <- cbind(KEY = seq_along(value), STATECENSUS = value)
    all <- merge(all, map)
    if (!is.na(year) && year != "NA") {
      if (isTRUE(startsWith(year, "mega"))) {
        suffix <- sub("^.*-([0-9]+)$", "\\1", year)
        if (suffix == "mega") {
        } else {
          n_years <- as.integer(suffix)
          all <- all[all$YEAR %in% seq(1994, 1994 + n_years - 1), ]
        }
        all$KEY <- (all$YEAR - min(all$YEAR)) * max(all$KEY) + all$KEY
      } else {
        year <- as.integer(year)
        all <- all[all$YEAR %in% year, ]
      }
    }
    all
  }
  all <- all[order(all$KEY), ]
  all
}

map <- function(all = NULL) {
  if (is.null(all)) all <- read_data()
  map <- unique(all[, c("KEY", "STATENAME", "STATECENSUS", "REGION", "DIVISION")])
  map <- map[order(map$KEY), ]
  map
}

prepare_for_dependent <- function(replication_number, shrinkage_shape, shrinkage_rate, grit_shape1, grit_shape2, years = "1994-2022", discard_proportion = 0.0) {
  if (grepl("[0-9]+_[0-9]+", replication_number)) {
    tmp <- strsplit(replication_number, "_", fixed = TRUE)[[1]]
    replication_number <- as.integer(tmp[2])
    k <- as.integer(tmp[1]) - 1L
  } else {
    replication_number <- as.integer(replication_number)
    k <- replication_number
  }
  anchor_anchor <- scan(file.path("anchors", paste0("regions", ".txt")), what = integer(), quiet = TRUE)
  anchor_shrinkage <- scan(file.path("shrinkages", paste0("borders", ".txt")), what = integer(), quiet = TRUE)
  anchor_shrinkage_reference <- 0
  all <- read_data()
  years <-if (years == "1994-2002") {
    1994:2002
  } else if (years == "1994-2011") {
    1994:2011
  } else if (years == "1994-2020") {
    1994:2020
  } else if (years == "1994-2020-2011") {
    c(1994:2020, 2020:2011)
  } else if (years == "1994-2020-2002") {
    c(1994:2020, 2020:2002)
  } else if (years == "1994-2020-1994") {
    c(1994:2020, 2020:1994)
  } else  {
    stop("Unsupported 'years' specification")
  }
  training_data_list <- vector("list", length(years))
  validation_data_list <- vector("list", length(years))
  for (i in seq_along(years)) {
    year <- years[i]
    all2 <- all[all$YEAR == year, ]
    map <- map(all2)
    n_items <- max(all2$KEY)
    training_which <- (seq_len(nrow(all2)) %% replication_number) != k
    training_data <- all2[training_which, , drop = FALSE]
    training_data <- training_data[sample(nrow(training_data), (1.0 - discard_proportion) * nrow(training_data)), , drop = FALSE]
    training_data <- training_data[order(training_data$KEY), , drop = FALSE]
    validation_data <- all2[!training_which, , drop = FALSE]
    if ((length(unique(training_data$KEY)) != n_items) || ((k != replication_number) && (length(unique(validation_data$KEY)) != n_items))) {
      stop("Bad validation and training split.")
    }
    training_data_list[[i]] <- training_data
    validation_data_list[[i]] <- validation_data
  }
  shrinkage_hyperparameters <-  list(reference = 1, shape = shrinkage_shape, rate = shrinkage_rate)
  grit_hyperparameters <- list(shape1 = grit_shape1, shape2 = grit_shape2)
  list(training_data_list = training_data_list,
       validation_data_list = if (k == replication_number) NULL else validation_data_list,
       shrinkage_hyperparameters = shrinkage_hyperparameters,
       grit_hyperparameters = grit_hyperparameters, anchor_anchor = anchor_anchor, anchor_shrinkage = anchor_shrinkage, anchor_shrinkage_reference = anchor_shrinkage_reference
      )
}

config <- function(year, partition_distribution, permutation, grit, anchor, shrinkage, replication_number, discard_proportion = 0.0) {
  if (grepl("[0-9]+_[0-9]+", replication_number)) {
    tmp <- strsplit(replication_number, "_", fixed = TRUE)[[1]]
    replication_number <- as.integer(tmp[2])
    k <- as.integer(tmp[1]) - 1L
  } else {
    replication_number <- as.integer(replication_number)
    k <- replication_number
  }
  # set.seed(replication_number)

  all <- read_data(year)
  map <- map(all)
  n_items <- max(all$KEY)

  training_which <- (seq_len(nrow(all)) %% replication_number) != k
  training_data <- all[training_which, , drop = FALSE]
  training_data <- training_data[sample(nrow(training_data), (1.0 - discard_proportion) * nrow(training_data)), , drop = FALSE]
  training_data <- training_data[order(training_data$KEY), , drop = FALSE]
  validation_data <- all[!training_which, , drop = FALSE]
  if ((length(unique(training_data$KEY)) != n_items) || ((k != replication_number) && (length(unique(validation_data$KEY)) != n_items))) {
    stop("Bad validation and training split.")
  }

  randomize_shrinkage <- FALSE
  shrinkage <- {
    shrinkage_hyperparameters = list(reference = 1L, shape = 5.0, rate = 1.0)
    shrinkage_file <- file.path("shrinkages", paste0(shrinkage, ".txt"))
    if (file.exists(shrinkage_file)) {
      scan(shrinkage_file, what = double(), quiet = TRUE)
    } else {
      if (!grepl("_", shrinkage, fixed = TRUE)) {
        as.double(shrinkage)
      } else {
        randomize_shrinkage <- TRUE
        tmp <- strsplit(shrinkage, "_", fixed = TRUE)[[1]]
        shrinkage_hyperparameters = list(reference = 1L, shape = as.double(tmp[1]), rate = as.double(tmp[2]))
        if (length(tmp) == 2) {
          shrinkage_hyperparameters$shape / shrinkage_hyperparameters$rate
        } else {
          shrinkage_hyperparameters$reference <- as.integer(tmp[3])
          scan(file.path("shrinkages", paste0(tmp[4], ".txt")), what = double(), quiet = TRUE)
        }
      }
    }
  }

  randomize_grit <- FALSE
  grit <- if (!grepl("_", grit, fixed = TRUE)) {
    grit_hyperparameters <- list(shape1 = 1.0, shape2 = 1.0)
    as.numeric(grit)
  } else {
    randomize_grit <- TRUE
    tmp <- as.numeric(strsplit(grit, "_", fixed = TRUE)[[1]])
    grit_hyperparameters <- list(shape1 = tmp[1], shape2 = tmp[2])
    grit_hyperparameters$shape1 / (grit_hyperparameters$shape1 + grit_hyperparameters$shape2)
  }

  anchor <- {
    option <- process_str(anchor)
    if (option[[1]]) {
      option[[2]]
    } else {
      if (anchor == "1") {
        rep(1, n_items)
      } else if (anchor == "n") {
        seq_len(n_items)
      } else {
        scan(file.path("anchors", paste0(anchor, ".txt")), what = integer(), quiet = TRUE)
      }
    }
  }

  randomize_permutation <- permutation == "random"
  permutation <- {
    if (permutation == "random") {
      map$KEY
    } else {
      p <- data.frame(PERMUTATION = seq_len(nrow(map)), STATENAME = readLines(paste0("permutations/", permutation, ".txt")))
      m <- merge(map, p, by = "STATENAME")
      m <- m[order(m$KEY), ]
      m$PERMUTATION
    }
  }

  is_crp_or_sp <- grepl("^(crp|sp_crp|sp_jl)[0-9]*[.]?[0-9]*$", partition_distribution)

  extract_concentration <- function() {
    concentration <- as.numeric(sub("^(crp|sp_crp|sp_jl)", "", partition_distribution))
    if (is.na(concentration)) concentration <- 1.0
    concentration
  }

  concentration <- 1
  discount <- 0.0
  baseline_crp <- CRPPartition(n_items, concentration = concentration, discount = discount)
  baseline_jl <- JensenLiuPartition(concentration = concentration, permutation = permutation)

  partition_distribution <- if (partition_distribution == "fixed") {
    FixedPartition(anchor)
  } else if (is_crp_or_sp && grepl("^crp", partition_distribution)) {
    CRPPartition(n_items, concentration = extract_concentration(), discount = discount)
  } else if (partition_distribution == "jl") {
    JensenLiuPartition(concentration = concentration, permutation = permutation)
  } else if (is_crp_or_sp && grepl("^sp_crp", partition_distribution)) {
    baseline_crp <- CRPPartition(n_items, concentration = extract_concentration(), discount = discount)
    ShrinkagePartition(anchor, shrinkage, permutation, grit, baseline_crp)
  } else if (is_crp_or_sp && grepl("^sp_jl", partition_distribution)) {
    baseline_jl <- JensenLiuPartition(concentration = extract_concentration(), permutation = permutation)
    ShrinkagePartition(anchor, shrinkage, permutation, grit, baseline_jl)
  } else if (partition_distribution == "lsp") {
    LocationScalePartition(anchor, shrinkage, concentration, permutation)
  } else if (partition_distribution == "cpp_crp_vi") {
    CenteredPartition(anchor, shrinkage, baseline_crp, useVI = TRUE)
  } else if (partition_distribution == "cpp_jl_vi") {
    CenteredPartition(anchor, shrinkage, baseline_jl, useVI = TRUE)
  } else if (partition_distribution == "cpp_crp_binder") {
    CenteredPartition(anchor, shrinkage, baseline_crp, useVI = FALSE)
  } else if (partition_distribution == "cpp_jl_binder") {
    CenteredPartition(anchor, shrinkage, baseline_jl, useVI = FALSE)
  } else {
    stop("Unrecognized partition distribution specification.")
  }

  list(
    partition_distribution = partition_distribution, randomize_permutation = randomize_permutation, randomize_shrinkage = randomize_shrinkage, randomize_grit = randomize_grit,
    grit_hyperparameters = grit_hyperparameters, shrinkage_hyperparameters = shrinkage_hyperparameters,
    anchor = anchor, n_items = n_items, map = map,
    training_data = training_data, validation_data = validation_data
  )
}

transform <- function(log_likelihood, n = n_reference, reference_log_likelihood = -55202.21507303197, n_reference = 139555) {
  1000000 * (( log_likelihood / n) - (reference_log_likelihood / n_reference))
}
