# Specific filter methods

# ------------------------------------------------------------------------------
# Area under the ROC curve

filter_roc_auc <-
  new_filter_method(
    name = "roc_auc",
    label = "Area under the ROC Curve",
    goal = "maximize",
    inputs = "quantitative",
    outputs = "qualitative",
    pkgs = "pROC"
  )

roc_wrapper <- function(x, y, ...) {
  cl <-
    rlang::call2(
      "roc_auc_vec",
      .ns = "yardstick",
      truth = quote(y),
      estimate = quote(x),
      ...
    )
  res <- try(rlang::eval_tidy(cl))
  if (inherits(res, "try-error")) {
    res <- NA_real_
  }
  res
}

#' @export
fit_xy.filter_method_roc_auc <- function(object, x, y, rename = FALSE, ...) {
  # check empty dots
  # case weights? event_level?
  x <- dplyr::as_tibble(x)
  if (is.vector(y)) {
    y <- dplyr::as_tibble(y)
  }
  validate_filter_data(object, x, y)

  p <- ncol(x)

  # Add wrapper using call and catch errors
  roc_scores <- purrr::map_dbl(x, ~ roc_wrapper(x = .x, y = y[[1]])) # TODO add in the ...
  roc_scores <- ifelse(roc_scores < 0.5, 1 - roc_scores, roc_scores)

  res <- new_filter_score(names(x), roc_scores, object, rename = rename,
                          num_pred = p)
  res
}

# ------------------------------------------------------------------------------
# Minimum Redundancy Maximal Relevancy Filter

filter_mrmr <-
  new_filter_method(
    name = "mrmr",
    label = "Minimum Redundancy Maximal Relevancy Filter",
    goal = "maximize",
    inputs = "all",
    outputs = "qualitative",
    pkgs = "praznik"
  )

#' @export
fit_xy.filter_method_mrmr <- function(object, x, y, rename = FALSE, ...) {
  x <- dplyr::as_tibble(x)
  if (is.vector(y)) { # TODO do these outside of these functions
    y <- dplyr::as_tibble(y)
  }
  validate_filter_data(object, x, y)
  # TODO convert ints in 'x' to doubles. See ?praznik::MRMR

  y <- y[[1]]
  p <- ncol(x)

  cl <-
    rlang::call2(
      "MRMR", .ns = "praznik",
      X = quote(x), Y = quote(y),
      k = p, ...
    )
  res <- try(rlang::eval_tidy(cl), silent = TRUE)

  if (inherits(res, "try-error")) {
    res$score <- rep(NA_real_, p)
    names(res) <- names(x)
  }

  res <- new_filter_score(names(res$score), unname(res$score), object,
                          rename = rename, num_pred = p)
  res
}

# ------------------------------------------------------------------------------
# Correlation filters


filter_corr <-
  new_filter_method(
    name = "corr",
    label = "Correlation Filter",
    goal = "maximize",
    inputs = "quantitative",
    outputs = "quantitative"
  )

#' @export
fit_xy.filter_method_corr <- function(object, x, y, rename = FALSE, ...) {
  x <- dplyr::as_tibble(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  p <- ncol(x)

  res <- cor(dplyr::bind_cols(y, x), use = "pairwise.complete.obs", ...)
  res <- abs(res[1, -1])

  res <- new_filter_score(names(x), res, object, rename = rename, num_pred = p)
  res
}

###

filter_rank_corr <-
  new_filter_method(
    name = "rank_corr",
    label = "Rank Correlation Filter",
    goal = "maximize",
    inputs = "quantitative",
    outputs = "quantitative"
  )

#' @export
fit_xy.filter_method_rank_corr <- function(object, x, y, rename = FALSE, ...) {
  fit_xy.filter_method_corr(object, x, y,  rename = rename, method = "spearman")
}

# ------------------------------------------------------------------------------
# Max difference in the outcome between groups of a factor predictor

filter_max_diff <-
  new_filter_method(
    name = "max_diff",
    label = "Maximum Group Difference",
    goal = "maximize",
    inputs = "qualitative",
    outputs = "all"
  )

comp_max_diff <- function(x, y) {
  dat <- data.frame(x = x, y = y)
  if (is.factor(y)) {
    fam <- "binomial"
  } else {
    fam <- "gaussian"
  }
  mod_fit <- try(stats::glm(y ~ x + 0, data = dat, family = fam), silent = TRUE)
  if (inherits(mod_fit, "try-error")) {
    res <- NA_real_
  } else {
    res <- coef(mod_fit)
    res <- max(res, na.rm = TRUE) - min(res, na.rm = TRUE)
  }
  res
}

#' @export
fit_xy.filter_method_max_diff <- function(object, x, y, rename = FALSE, ...) {
  x <- dplyr::as_tibble(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  p <- ncol(x)

  res <- purrr::map_dbl(x, ~ comp_max_diff(.x, y = y[[1]]))
  res <- new_filter_score(names(x), res, object, rename = rename, num_pred = p)
  res
}

# ------------------------------------------------------------------------------
# Variable importance

filter_rf_imp <-
  new_filter_method(
    name = "rf_imp",
    label = "Random Forest Variable Importance",
    goal = "maximize",
    inputs = "all",
    outputs = "all",
    pkgs = "ranger"
  )

#' @export
fit_xy.filter_method_rf_imp <- function(object, x, y, rename = FALSE,
                                        seed = sample.int(1000, 1), ...) {
  x <- dplyr::as_tibble(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  y <- y[[1]]
  p <- ncol(x)

  # check dots for 'importance'; when not set or use 'none'
  cl <- rlang::call2("ranger", .ns = "ranger", x = quote(x), y = quote(y),
                     importance = "impurity_corrected", ...)
  set.seed(seed) # TODO use withr::with_seed?
  res <- try(rlang::eval_tidy(cl), silent = TRUE)

  if (inherits(res, "try-error")) {
    res <- list(variable.importance = rep(NA_real_, p))
    names(res$variable.importance) <- names(x)
  }

  res <-
    new_filter_score(
      names(res$variable.importance),
      unname(res$variable.importance),
      object,
      rename = rename,
      num_pred = p
    )
  res
}


# ------------------------------------------------------------------------------
# Information


filter_info_gain <-
  new_filter_method(
    name = "info_gain",
    label = "Information Gain",
    goal = "maximize",
    inputs = "all",
    outputs = "qualitative",
    pkgs = "FSelectorRcpp"
  )


#' @export
fit_xy.filter_method_info_gain <- function(object, x, y, rename = FALSE, ...) {
  x <- as.data.frame(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  y <- y[[1]]
  p <- ncol(x)

  cl <- rlang::call2("information_gain", .ns = "FSelectorRcpp",
                     x = quote(x), y = quote(y), type = "infogain", ...)

  res <- try(rlang::eval_tidy(cl), silent = TRUE)

  if (inherits(res, "try-error")) {
    res <- dplyr::tibble(variable = names(x), score = rep(NA_real_, p))
  } else {
    res <- setNames(res, c("variable", "score"))
  }

  res <-
    new_filter_score(
      res$variable,
      res$score,
      object,
      rename = rename,
      num_pred = p
    )
  res
}

###

filter_info_gain_ratio <-
  new_filter_method(
    name = "info_gain_ratio",
    label = "Information Gain Ratio",
    goal = "maximize",
    inputs = "all",
    outputs = "qualitative",
    pkgs = "FSelectorRcpp"
  )

#' @export
fit_xy.filter_method_info_gain_ratio <- function(object, x, y, rename = FALSE, ...) {
  x <- as.data.frame(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  y <- y[[1]]
  p <- ncol(x)

  cl <- rlang::call2("information_gain", .ns = "FSelectorRcpp",
                     x = quote(x), y = quote(y), type = "gainratio", ...)

  res <- try(rlang::eval_tidy(cl), silent = TRUE)

  if (inherits(res, "try-error")) {
    res <- dplyr::tibble(variable = names(x), score = rep(NA_real_, p))
  } else {
    res <- setNames(res, c("variable", "score"))
  }

  res <-
    new_filter_score(
      res$variable,
      res$score,
      object,
      rename = rename,
      num_pred = p
    )
  res
}

###

filter_mic <-
  new_filter_method(
    name = "mic",
    label = "Maximal Information Coefficient",
    goal = "maximize",
    inputs = "quantitative",
    outputs = "quantitative",
    pkgs = "minerva"
  )

#' @export
fit_xy.filter_method_mic <- function(object, x, y, rename = FALSE, ...) {
  x <- dplyr::as_tibble(x)
  y <- dplyr::as_tibble(y)
  validate_filter_data(object, x, y)

  x <- as.matrix(x)
  y <- y[[1]]
  p <- ncol(x)

  cl <- rlang::call2("mine", .ns = "minerva",
                     x = quote(x), y = quote(y),
                     use = "pairwise.complete.obs", ...)

  res <- try(rlang::eval_tidy(cl), silent = TRUE)

  if (inherits(res, "try-error")) {
    res <- rep(NA_real_, p)
  } else {
    res <- res$MIC[,1]
  }

  res <-
    new_filter_score(
      colnames(x),
      res,
      object,
      rename = rename,
      num_pred = p
    )
  res
}

