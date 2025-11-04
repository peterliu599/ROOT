# ==========================================================
#' Check for missing values in training data
# ==========================================================
#'
#' Ensures there are no NA values in any of the relevant columns of the dataset.
#'
#' @param data A data frame.
#' @param cols Character vector of column names to check.
#' @return Invisibly returns TRUE if no NA found; otherwise throws an error.
#' @keywords internal
check_no_na <- function(data, cols) {
  data_name <- deparse(substitute(data))
  for (col in cols) {
    if (anyNA(data[[col]])) {
      stop(sprintf("Data `%s` column '%s' contains missing values. Please handle NA before training.", data_name, col),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ==========================================================
#' Train nuisance models for weighting
# ==========================================================
#'
#' Fits models to estimate sampling and treatment propensities on training data.
#' Specifically, it computes:
#' \itemize{
#'   \item \code{pi}: The prevalence of sample inclusion (estimated as mean of `S`).
#'   \item \code{pi_m}: A logistic regression model for \eqn{P(S=1 \mid X)} using all covariates.
#'   \item \code{e_m}: A logistic regression model for \eqn{P(Tr=1 \mid X, S=1)}, fit only on the subset where `S==1`.
#' }
#'
#' @param training_data A data frame containing the training dataset.
#' @param outcome Name of the outcome column (typically observed outcome, e.g. "Yobs").
#' @param treatment Name of the treatment indicator column (e.g. "Tr").
#' @param sample Name of the sample inclusion indicator column (e.g. "S").
#' @return A list with components:
#' \item{pi}{Numeric scalar giving the overall sample inclusion rate \eqn{P(S=1)} in the training data.}
#' \item{pi_m}{A fitted `glm` model (binomial family) for \eqn{P(S=1 \mid X)}.}
#' \item{e_m}{A fitted `glm` model (binomial family) for \eqn{P(Tr=1 \mid X, S=1)}.}
#' @details This function uses simple logistic regression (`glm` with logit link) to estimate the necessary nuisance parameters for weighting.
#'   It requires that both `S` and `Tr` have variation (both 0 and 1 must be present in training data); if not, the fitting is not possible and an error is raised.
#'   All covariates other than the specified outcome, treatment, and sample columns are used as predictors.
train <- function(training_data, outcome, treatment, sample) {
  # Input validation
  if (!is.data.frame(training_data)) {
    stop("`training_data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(training_data)) {
      stop(sprintf("Column '%s' not found in training_data.", col), call. = FALSE)
    }
  }
  check_no_na(training_data, colnames(training_data))
  # Ensure binary columns are numeric 0/1
  S_vec <- training_data[[sample]]
  Tr_vec <- training_data[[treatment]]
  if (!is.numeric(S_vec) || !all(S_vec %in% c(0, 1))) {
    stop("Sample indicator column must be numeric 0/1.", call. = FALSE)
  }
  if (!is.numeric(Tr_vec) || !all(Tr_vec %in% c(0, 1))) {
    stop("Treatment indicator column must be numeric 0/1.", call. = FALSE)
  }
  # Require both classes present
  if (all(S_vec == 0) || all(S_vec == 1)) {
    stop("Cannot train model: sample indicator column `", sample, "` has no variation (all values are the same).", call. = FALSE)
  }
  if (all(Tr_vec[S_vec == 1] == 0) || all(Tr_vec[S_vec == 1] == 1)) {
    stop("Cannot train model: treatment column has no variation among S==1 observations.", call. = FALSE)
  }

  # Partition covariates (X) versus specified columns
  covariate_cols <- setdiff(names(training_data), c(outcome, treatment, sample))
  X_df <- training_data[, covariate_cols, drop = FALSE]

  # Compute sample prevalence
  pi <- mean(S_vec)

  # Fit logistic model for P(S=1 | X)
  data_pi <- data.frame(S = S_vec, X_df)
  pi_m <- stats::glm(formula = S ~ ., data = data_pi, family = stats::binomial())

  # Fit logistic model for P(Tr=1 | X, S=1) using only S==1 subset
  data_e <- training_data[training_data[[sample]] == 1, c(treatment, covariate_cols), drop = FALSE]
  e_m <- stats::glm(formula = stats::as.formula(paste(treatment, "~ .")), data = data_e, family = stats::binomial())

  return(list(pi = pi, pi_m = pi_m, e_m = e_m))
}

#' Compute pseudo-outcome components (a, b) and their product (v)
#'
#' Using the outputs of the nuisance models, computes intermediate values for the
#' treatment effect estimation via inverse probability weighting (IPW) for the ATT in sample.
#'
#' Specifically, it computes:
#' \itemize{
#'   \item \code{a}: IPW-adjusted outcome difference, \eqn{a_i = S_i \left(\frac{Tr_i Y_i}{p_{t1|x,i}} - \frac{(1-Tr_i) Y_i}{1 - p_{t1|x,i}}\right)}.
#'   \item \code{b}: Overlap weight factor, \eqn{b_i = \frac{1}{\ell(X_i)}}, where \eqn{\ell(X) = \frac{P(S=1|X)/\pi}{P(S=0|X)/(1-\pi)}}.
#'   \item \code{v}: The pseudo-outcome, defined as \eqn{v_i = a_i \times b_i}.
#' }
#'
#' @param testing_data A data frame of test data (or evaluation data) containing at least the columns for outcome, treatment, and sample indicators.
#' @param outcome Name of the outcome column in `testing_data`.
#' @param treatment Name of the treatment column in `testing_data` (0/1).
#' @param sample Name of the sample indicator column in `testing_data` (0/1).
#' @param pi Numeric scalar, the estimated \eqn{P(S=1)} (prevalence) from the training data.
#' @param pi_m A fitted model (e.g., `glm`) for \eqn{P(S=1 \mid X)}; typically from `train()`.
#' @param e_m A fitted model (`glm`) for \eqn{P(Tr=1 \mid X, S=1)}; typically from `train()`.
#' @return A list with numeric vectors:
#' \item{v}{Pseudo-outcome values for each observation (numeric vector length = nrow(testing_data)).}
#' \item{a}{Intermediate "IPW-adjusted outcome" values (same length as v).}
#' \item{b}{Overlap weight factors (same length as v).}
#' @details The predicted probabilities from `pi_m` and `e_m` are constrained to `[1e-8, 1-1e-8]` to avoid instability (extremely small or large probabilities are clamped).
#'   If the provided `pi` is 0 or 1 (indicating no variation in sample inclusion in training), the computation is undefined and an error will be thrown.
#'   Ensure that `pi_m` and `e_m` correspond to models trained on compatible data (same covariates) for accurate predictions.
#' @seealso \code{\link{train}} for obtaining `pi`, `pi_m`, and `e_m`; \code{\link{estimate_dml}} for cross-fitted estimation.
estimate <- function(testing_data, outcome, treatment, sample, pi, pi_m, e_m) {
  # Input validation
  if (!is.data.frame(testing_data)) {
    stop("`testing_data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(testing_data)) {
      stop(sprintf("Column '%s' not found in testing_data.", col), call. = FALSE)
    }
  }
  check_no_na(testing_data, colnames(testing_data))
  if (!is.numeric(pi) || length(pi) != 1) {
    stop("`pi` must be a numeric scalar.", call. = FALSE)
  }
  if (pi <= 0 || pi >= 1) {
    stop("Invalid `pi`: P(S=1) must be between 0 and 1 (exclusive).", call. = FALSE)
  }
  if (!inherits(pi_m, "glm") || !inherits(e_m, "glm")) {
    stop("`pi_m` and `e_m` must be model objects (e.g. from glm) for prediction.", call. = FALSE)
  }

  # Extract columns from testing_data
  S_vec  <- testing_data[[sample]]
  Y_vec  <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  # Data frame of covariates (all other columns)
  covariate_cols <- setdiff(names(testing_data), c(outcome, treatment, sample))
  X_df <- testing_data[, covariate_cols, drop = FALSE]

  # Predict P(S=1 | X) and compute overlap factor lX
  p_s1x <- stats::predict(pi_m, newdata = X_df, type = "response")
  # Clamp predicted probabilities to avoid 0 or 1
  p_s1x <- pmin(pmax(as.numeric(p_s1x), 1e-8), 1 - 1e-8)
  lX <- (p_s1x / pi) / ((1 - p_s1x) / (1 - pi))

  # Predict P(Tr=1 | X, S=1) for all observations (treat those with S=0 similarly)
  p_t1x <- stats::predict(e_m, newdata = X_df, type = "response")
  p_t1x <- pmin(pmax(as.numeric(p_t1x), 1e-8), 1 - 1e-8)

  # Compute IPW components
  a_val <- S_vec * (Tr_vec * Y_vec / p_t1x - (1 - Tr_vec) * Y_vec / (1 - p_t1x))
  b_val <- 1 / lX
  v_val <- a_val * b_val

  return(list(v = as.numeric(v_val), a = as.numeric(a_val), b = as.numeric(b_val)))
}

#' Stratified K-fold index generator
#'
#' Splits indices into K folds while preserving the class distribution of a binary factor.
#' This mimics scikit-learn's `StratifiedKFold`, ensuring each fold has a representative ratio of the two classes in `S`.
#'
#' @param S A vector or factor indicating class membership (typically 0/1 or two-class factor) for stratification.
#' @param K Integer number of folds (K >= 2). If K is larger than the number of observations, it will be reduced to that number.
#' @return A list of length K, where each element is an integer vector of row indices assigned to that fold.
#'   The union of all folds equals \code{1:length(S)}, and folds are roughly equal in size.
#' @details The function deterministically allocates indices to folds by class. For each class in `S`, indices are cyclically assigned to folds to balance counts.
#'   If `K == 1`, a single fold containing all indices is returned (though typically K should be >= 2 for cross-validation).
stratified_kfold <- function(S, K = 5) {
  # Input validation
  if (is.data.frame(S) || is.matrix(S)) {
    stop("`S` should be a vector or factor, not a data frame or matrix.", call. = FALSE)
  }
  if (!is.factor(S)) {
    # Coerce S to factor (works for numeric 0/1 as well)
    S <- as.factor(S)
  }
  if (!is.numeric(K) || length(K) != 1 || K < 1) {
    stop("`K` must be a positive integer.", call. = FALSE)
  }
  K <- as.integer(K)
  if (K > length(S)) {
    warning("Requested K (", K, ") is greater than number of observations; using K = ", length(S), ".", call. = FALSE)
    K <- length(S)
  }
  if (K < 1) {
    K <- 1  # edge-case safeguard
  }

  # Split indices by class
  idx_by_class <- split(seq_along(S), S)
  # Distribute each class's indices across K folds in round-robin fashion
  parts_by_class <- lapply(idx_by_class, function(idx) {
    split(idx, rep_len(seq_len(K), length(idx)))
  })

  # Initialize list of folds
  folds <- vector("list", K)
  for (k in seq_len(K)) {
    # Combine indices from each class for fold k
    fold_k_indices <- unlist(lapply(parts_by_class, `[[`, k), use.names = FALSE)
    folds[[k]] <- sort(fold_k_indices)
  }
  return(folds)
}

#' Cross-fitted estimation of pseudo-outcomes (Double ML)
#'
#' Trains nuisance models on each training fold and computes pseudo-outcomes on
#' the corresponding test fold, then aggregates results. Returns only what is
#' needed downstream: the pseudo-outcome table and the aligned evaluation data.
#'
#' @param data A data frame containing at least the outcome, treatment, and sample indicator columns.
#' @param outcome Name of the outcome column.
#' @param treatment Name of the treatment column (0/1).
#' @param sample Name of the sample indicator column (0/1).
#' @param crossfit Integer number of folds for cross-fitting (>= 2).
#' @return A list with:
#'   \item{df_v}{Data frame with one row per kept observation (indexed by `primary_index`),
#'               containing: `te` (pseudo-outcome v), `a`, `b`, and squared deviations
#'               `te_sq`, `a_sq`. Only S==1 rows with finite values are kept.}
#'   \item{data2}{Subset of original `data` corresponding to `df_v$primary_index`.}
#' @note Rows with infinite or undefined weights (e.g., where the predicted propensity scores were 0 or 1) are removed from `df_v` (and the corresponding rows in `data2`). The `primary_index` in `df_v` corresponds to the row index in the original `data`.
#'   Squared deviation columns (`te_sq`, `a_sq`) are centered around the mean of `te` and `a` for the `S==1` group.
#' @examples
#'  sim<-get_data(n=2000,seed=599)
#'  dml<-estimate_dml(sim$data,outcome="Yobs",treatment="Tr",sample="S",crossfit= 5)
#' @export
estimate_dml <- function(data, outcome, treatment, sample, crossfit = 5) {
  # Input validation
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(data)) stop(sprintf("Column '%s' not found in data.", col), call. = FALSE)
  }
  check_no_na(data, colnames(data))
  if (!is.numeric(crossfit) || length(crossfit) != 1 || crossfit < 2) {
    stop("`crossfit` must be an integer >= 2.", call. = FALSE)
  }
  crossfit <- as.integer(crossfit)
  if (all(data[[sample]] == 0) || all(data[[sample]] == 1)) {
    stop("Sample indicator `", sample, "` has no variation (all 0 or all 1). Cannot perform cross-fitting.", call. = FALSE)
  }

  n <- nrow(data)
  folds <- stratified_kfold(data[[sample]], K = crossfit)

  df_list <- vector("list", length(folds))

  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)

    training_data <- data[train_idx, , drop = FALSE]
    testing_data  <- data[test_idx,  , drop = FALSE]

    # Train nuisances on training fold
    trn <- train(training_data, outcome = outcome, treatment = treatment, sample = sample)

    # Pseudo-outcomes on test fold
    est <- estimate(testing_data,
                    outcome   = outcome,
                    treatment = treatment,
                    sample    = sample,
                    pi        = trn$pi,
                    pi_m      = trn$pi_m,
                    e_m       = trn$e_m)

    df_list[[i]] <- data.frame(
      te            = as.numeric(est$v),
      primary_index = test_idx,
      a             = as.numeric(est$a),
      b             = as.numeric(est$b),
      row.names     = NULL
    )
  }

  # Combine folds
  df_v_all <- do.call(rbind, df_list)

  # Handle Inf/NA and drop bad rows
  is_inf <- sapply(df_v_all, is.infinite)
  if (any(is_inf)) df_v_all[is_inf] <- NA
  df_v_all <- stats::na.omit(df_v_all)

  # Means of te and a over S==1 group
  s_idx1 <- which(data[[sample]] == 1)
  te_mean_s1 <- mean(df_v_all$te[df_v_all$primary_index %in% s_idx1])
  a_mean_s1  <- mean(df_v_all$a[df_v_all$primary_index %in% s_idx1])

  # Add squared deviation columns
  df_v_all$te_sq <- (df_v_all$te - te_mean_s1)^2
  df_v_all$a_sq  <- (df_v_all$a - a_mean_s1)^2

  # Aggregate (each index appears once with proper K-folds, but safe to aggregate)
  df_v_grp <- stats::aggregate(df_v_all[, c("te", "a", "b", "te_sq", "a_sq")],
                        by = list(primary_index = df_v_all$primary_index),
                        FUN = mean)

  # Keep only S==1 indices
  df_v_grp <- df_v_grp[df_v_grp$primary_index %in% s_idx1, , drop = FALSE]

  # Align original data rows
  data_s1 <- data[df_v_grp$primary_index, , drop = FALSE]

  list(df_v = df_v_grp, data2 = data_s1)
}


#' Fit a shallow decision tree to characterize learned weights `w`
#'
#' Trains a classification tree on the covariates `X` to predict the binary membership `w`.
#' This provides an interpretable summary of how the weighted subgroup can be distinguished by `X`.
#'
#' @param X A data frame of covariates (features).
#' @param w A binary vector (0/1 or a factor with two levels) indicating class membership for each observation (e.g., whether an observation is in the selected subgroup).
#' @param max_depth Integer, the maximum tree depth (default 3).
#' @return An `rpart` object representing the fitted decision tree.
#' @details The tree is grown using the Gini index (classification) and is not pruned (complexity parameter `cp = 0`), relying solely on `max_depth` to control complexity. This mirrors the default behavior of scikit-learn's `DecisionTreeClassifier(max_depth=...)`.
#'   If `w` is not already a factor, it will be converted internally. The tree's rules can be interpreted to understand which covariates (and what splits) best separate the two classes defined by `w`.
characterize_tree <- function(X, w, max_depth = 3) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  check_no_na(X, colnames(X))
  if (length(w) != nrow(X)) {
    stop("Length of `w` must equal the number of rows in `X`.", call. = FALSE)
  }
  # Coerce w to factor and check it has two levels
  w_factor <- as.factor(w)
  if (nlevels(w_factor) != 2) {
    stop("`w` must have exactly two classes (binary).", call. = FALSE)
  }

  # Prepare data for rpart
  df <- data.frame(w = w_factor, X)
  # Fit classification tree
  fit <- rpart::rpart(
    w ~ .,
    data   = df,
    method = "class",
    parms  = list(split = "gini"),
    control = rpart::rpart.control(maxdepth = max_depth, cp = 0.0, minsplit = 2, minbucket = 1),
    model = TRUE
  )
  return(fit)
}

#' Compute the midpoint of a numeric vector
#'
#' Calculates the midpoint defined as \eqn{( \max(X) + \min(X) ) / 2}, ignoring any NA values.
#'
#' @param X A numeric vector.
#' @return A numeric scalar giving the midpoint of the finite values in `X`. If `X` is empty or has no finite values, `NA` is returned.
midpoint <- function(X) {
  # Input validation
  if (!is.numeric(X)) {
    stop("`X` must be a numeric vector.", call. = FALSE)
  }
  if (length(X) == 0) {
    warning("Empty vector provided to midpoint(). Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  check_no_na(X, colnames(X))
  # Remove NAs for calculation
  rng <- range(X, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    # All values are NA or infinite
    warning("No finite values in `X`. Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  return((rng[1] + rng[2]) / 2)
}

#' Randomly choose a split feature based on provided probabilities
#'
#' Given a probability distribution over features (and possibly a "leaf" option), selects one feature at random according to those probabilities.
#'
#' @param split_feature A named numeric vector of feature selection probabilities. Names should correspond to feature IDs (and may include a special "leaf" entry).
#' @param depth Current tree depth (an integer, used for parity with Python implementation but not affecting probabilities in this implementation).
#' @return A single feature name (or "leaf") chosen randomly according to the provided probability weights.
#' @note The factor \eqn{2^{(0 * depth / 4)}} present in the code is effectively 1 (no effect on the first element's weight) and is included only for parity with an equivalent Python implementation. All probabilities are normalized to sum to 1 before sampling.
choose_feature <- function(split_feature, depth) {
  # Input validation
  if (!is.numeric(split_feature) || length(split_feature) == 0) {
    stop("`split_feature` must be a non-empty numeric vector of probabilities.", call. = FALSE)
  }
  if (is.null(names(split_feature))) {
    stop("`split_feature` must have names for each feature (use 'leaf' for the special leaf option if applicable).", call. = FALSE)
  }
  if (!is.numeric(depth) || length(depth) != 1) {
    stop("`depth` must be a numeric scalar (current tree depth).", call. = FALSE)
  }

  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf("All probabilities must be >= 0. Negative entries found in: %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }

  split_prob <- as.numeric(split_feature)
  # Adjust the first element probability by the (inactive) depth factor (kept for parity)
  if (length(split_prob) > 0) {
    split_prob[1] <- split_prob[1] * (2^(0 * depth / 4))
  }
  # Normalize to a probability distribution
  total <- sum(split_prob)
  if (!is.finite(total)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }
  split_prob <- split_prob / total

  # Perform random sampling
  choices <- names(split_feature)
  chosen <- sample(choices, size = 1, replace = TRUE, prob = split_prob)
  return(chosen)
}

#' Evaluate splitting objective loss for a given weight assignment
#'
#' Computes the approximate standard error (loss) of the treatment effect estimate under a hypothetical assignment of weights `w` for specified rows.
#'
#' @param val Numeric scalar (must be 0 or 1). The weight value to assign (0 = exclude, 1 = include).
#' @param indices Indices or row names in `D` for which the weight should be set to `val`. Can be a numeric vector of row positions or a character vector of row names.
#' @param D A data frame containing at least columns `vsq` (squared pseudo-outcome) and `w` (current weights).
#' @return Numeric value representing the loss, defined as \eqn{\sqrt{\sum_i vsq_i * w_i \,/\, (\sum_i w_i)^2}}. Returns `Inf` if the denominator is 0 or if the result is not a number.
#' @note This function mimics the behavior of numpy's `nan_to_num(..., nan=Inf)` by returning `Inf` when the computation is undefined (e.g., no weights selected). It is used internally to decide whether a proposed split improves the objective.
#' @examples
#' # Create a small example data frame D
#' D <- data.frame(vsq = c(0.5, 1.0, 0.2), w = c(1, 1, 1))
#' loss(0, indices = 2, D)  # Set w=0 for row 2 and compute loss
#' loss(1, indices = c(2,3), D)  # Set w=1 for rows 2 and 3 and compute loss
#' @export
loss <- function(val, indices, D) {
  # Input validation
  if (!is.numeric(val) || length(val) != 1 || !(val %in% c(0, 1))) {
    stop("`val` must be a numeric scalar equal to 0 or 1.", call. = FALSE)
  }
  if (!is.data.frame(D) || !all(c("vsq", "w") %in% names(D))) {
    stop("`D` must be a data frame with numeric columns 'vsq' and 'w'.", call. = FALSE)
  }
  # Determine rows to update
  rows_to_update <- integer(0)
  if (length(indices) > 0) {
    if (is.numeric(indices)) {
      rows_to_update <- as.integer(indices)
    } else {
      # assume character or factor -> treat as row names
      rows_to_update <- which(rownames(D) %in% as.character(indices))
    }
  }
  # Copy D to avoid modifying original
  D_temp <- D
  if (length(rows_to_update) > 0) {
    D_temp[rows_to_update, "w"] <- val
  }
  # Compute numerator and denominator for loss
  num <- sum(D_temp$vsq * D_temp$w, na.rm = TRUE)
  den <- (sum(D_temp$w, na.rm = TRUE))^2
  loss_val <- sqrt(num / den)
  if (!is.finite(loss_val) || is.nan(loss_val)) {
    loss_val <- Inf
  }
  return(loss_val)
}



#' Reduce a feature's selection weight by half and renormalize
#'
#' Lowers the probability weight of a given feature by 50%, and then re-normalizes the entire probability vector.
#'
#' @param fj A feature name (character string) present in the names of `split_feature`.
#' @param split_feature A named numeric vector of probabilities for features (as used in splitting).
#' @return A numeric vector of the same length as `split_feature`, giving the updated probabilities that sum to 1.
#' @details This is typically used when a particular feature split was rejected; the feature's probability is halved to reduce its chance of being chosen again immediately, encouraging exploration of other features. If `fj` is "leaf", its weight is also halved similarly.
reduce_weight <- function(fj, split_feature) {
  # Input validation
  if (!is.character(fj) || length(fj) != 1) {
    stop("`fj` must be a single feature name (string).", call. = FALSE)
  }
  if (!is.numeric(split_feature) || is.null(names(split_feature))) {
    stop("`split_feature` must be a named numeric vector of probabilities.", call. = FALSE)
  }
  if (!(fj %in% names(split_feature))) {
    stop("Feature '", fj, "' not found in names(split_feature).", call. = FALSE)
  }
  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf("All probabilities must be >= 0. Negative entries found in: %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }

  total_before <- sum(split_feature)
  if (!is.finite(total_before)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }

  # Reduce the specified feature's weight by half
  split_feature[fj] <- split_feature[fj] / 2
  # Renormalize to sum to 1
  split_feature <- split_feature / sum(split_feature)
  return(split_feature)
}

#' Recursive split builder for weighted tree (internal function)
#'
#' Recursively builds a weighted decision tree to optimize an objective function, using exploration/exploitation trade-off. This function is internal and is used by `tree_opt()` to construct a single tree.
#'
#' @param split_feature Named numeric vector of feature selection probabilities (should include a "leaf" option).
#' @param X A data frame of current observations (must include at least the feature `fj` chosen for splitting if applicable; may also include a working copy of weights `w`).
#' @param D A data frame representing the global state (must include columns `w` for weights and `vsq` for squared outcomes, with row names aligning to observations).
#' @param parent_loss Numeric, the loss value of the parent node (used to decide if a split improves the objective).
#' @param depth Integer, current depth of the node in the tree.
#' @param explore_proba Numeric in `[0,1]`, probability of exploring the suboptimal split at a leaf (randomly flipping the chosen weight).
#' @param choose_feature Function to choose the next feature to split (default uses `choose_feature`).
#' @param loss_fn Function to compute loss given a tentative weight assignment (default uses `loss`).
#' @param reduce_weight_fn Function to adjust feature weights when a split is rejected (default uses `reduce_weight`).
#' @param max_depth Integer maximum depth of the tree. If the current depth equals `max_depth`, the node is made a leaf.
#' @param min_leaf_n Integer minimum number of observations required to attempt a split. If `X` has <= `min_leaf_n` rows, the node becomes a leaf.
#' @param log_fn Function for logging debug messages. By default, a no-op. If provided (e.g., from `tree_opt(verbose=TRUE)`), it will receive formatted strings to output.
#' @return A list representing the tree/subtree. Each node includes:
#'   \item{node}{Feature name used for split at this node (or "leaf").}
#'   \item{split}{Numeric value of the splitting threshold (midpoint) if node is a split.}
#'   \item{left_tree}{Left subtree (list) if node is a split.}
#'   \item{right_tree}{Right subtree (list) if node is a split.}
#'   \item{w}{If node is a leaf, the weight (0 or 1) assigned to all observations in this node.}
#'   \item{local objective}{Numeric, the loss objective value at this node after assignment.}
#'   \item{depth}{Integer depth of this node.}
#'   \item{D}{Data frame of global state after this node's assignments (for internal use).}
#'   \item{leaf_reason}{(Leaf nodes only) Text reason for why this node was made a leaf (e.g., "max-depth", "min-leaf", "feature==leaf", etc.).}
#'   \item{feature}{(Leaf nodes only) The feature that would have been split (or NA).}
#'   \item{cut}{(Leaf nodes only) The cut value (or NA) that would have been used.}
#' @details This function works as follows:
#'   \enumerate{
#'     \item Selects a feature `fj` to split using `choose_feature` and the current probability vector `split_feature`. If `fj == "leaf"`, no further splitting is done and the node becomes a leaf.
#'     \item If the maximum depth or minimum node size criteria are met, the node becomes a leaf.
#'     \item Otherwise, it computes a split threshold (`cj`) as the midpoint of feature `fj` in the current data `X`, and partitions `X` (and the corresponding indices in `D`) into left (`X_left`) and right (`X_right`) sets.
#'     \item If either side is empty, the node becomes a leaf (no split possible).
#'     \item It then evaluates the loss for assigning all left observations weight 0 vs 1, and all right observations weight 0 vs 1, to determine optimal assignments (`w_left` and `w_right` for each side).
#'     \item If the combined loss `new_loss` is less than or equal to the `parent_loss`, the split is accepted: weights in `D` are updated for left and right, and the function recurses on each side (randomizing the order of recursion to break symmetry).
#'     \item If the split does not improve the loss, it is rejected: the selection probability of `fj` is reduced (via `reduce_weight_fn`), and `split_node` is called again on the same data without increasing depth.
#'   }
#'   Leaf node assignment: when a node is designated as a leaf, it compares the loss of assigning all observations in that node `w=0` versus `w=1` and picks the assignment that yields lower loss (exploitation). With probability `explore_proba`, it may instead pick the opposite assignment to encourage exploration. This randomness can be controlled via the `explore_proba` parameter.
#' @note This function is typically not called directly by the user. It assumes that `D` and `X` are kept in sync and that `D$w` is updated globally as splits are made. The `log_fn` parameter allows optional logging of the splitting process for debugging or insight when `verbose=TRUE` is passed from higher-level functions.
split_node <- function(split_feature, X, D, parent_loss, depth,
                       explore_proba = 0.05,
                       choose_feature = choose_feature,
                       loss_fn = loss,
                       reduce_weight_fn = reduce_weight,
                       max_depth = 8,
                       min_leaf_n = 5,
                       log_fn = function(...) {}) {
  # Input validation
  if (!is.numeric(parent_loss) || length(parent_loss) != 1) {
    stop("`parent_loss` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(depth) || length(depth) != 1) {
    stop("`depth` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(explore_proba) || explore_proba < 0 || explore_proba > 1) {
    stop("`explore_proba` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.data.frame(X) || !is.data.frame(D)) {
    stop("`X` and `D` must be data frames.", call. = FALSE)
  }
  check_no_na(X, colnames(X))
  check_no_na(D, colnames(D))
  if (!"w" %in% names(D) || !"vsq" %in% names(D)) {
    stop("`D` must contain columns 'w' and 'vsq'.", call. = FALSE)
  }
  if (!is.numeric(split_feature) || is.null(names(split_feature)) || length(split_feature) == 0) {
    stop("`split_feature` must be a named numeric vector (including 'leaf').", call. = FALSE)
  }
  if (!"leaf" %in% names(split_feature)) {
    stop("`split_feature` must include an entry named 'leaf'.", call. = FALSE)
  }

  # Define a helper logging function
  .log <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    log_fn(msg)
  }
  # Define a nearly-less-or-equal for floating comparisons (tolerance for improvement)
  nearly_leq <- function(a, b, tol_abs = 1e-12, tol_rel = 1e-12) {
    a <= b + tol_abs + tol_rel * max(1, abs(b))
  }
  # Define how to create a leaf node result
  make_leaf <- function(X_sub, D_sub, depth_cur, reason, fj = NA_character_, cj = NA_real_) {
    # Compute losses for assigning all in this node w=0 vs w=1
    losses <- c(loss_fn(0, rownames(X_sub), D_sub), loss_fn(1, rownames(X_sub), D_sub))
    w_exploit <- which.min(losses) - 1  # optimal assignment (0 or 1 that minimizes loss)
    w_explore <- stats::rbinom(1, 1, 0.5)  # random 0 or 1
    explore_flip <- stats::rbinom(1, 1, explore_proba)
    # Choose w: with probability explore_proba, use the random choice, otherwise use optimal
    final_w <- if (explore_flip == 1) w_explore else w_exploit


    .log("[depth=%d] LEAF (%s): feature=%s, cut=%s, n=%d, losses={%.4f, %.4f}, w=%d",
         depth_cur, reason, as.character(fj),
         ifelse(is.na(cj), "NA", sprintf("%.4f", cj)),
         nrow(X_sub), losses[1], losses[2], final_w)


    # Assign weights in D for this node's observations
    if (nrow(X_sub) > 0) {
      D_sub[rownames(X_sub), "w"] <- final_w
      X_sub$w <- final_w
    }
    return(list(node = "leaf", w = final_w,
                `local objective` = min(losses), depth = depth_cur, D = D_sub,
                leaf_reason = reason, feature = fj, cut = cj))
  }

  # Stopping rule 1: max depth reached
  if (depth >= max_depth) {
    return(make_leaf(X, D, depth, reason = "max-depth", fj = NA_character_))
  }
  # Stopping rule 2: too few observations to split
  if (nrow(X) <= min_leaf_n) {
    return(make_leaf(X, D, depth, reason = "min-leaf", fj = NA_character_))
  }

  # Choose a feature to split (could be "leaf")
  fj <- choose_feature(split_feature, depth)

  # If the chosen feature is "leaf", make this node a leaf (stop splitting)
  if (identical(fj, "leaf")) {
    return(make_leaf(X, D, depth, reason = "feature==leaf", fj = fj))
  }

  # Determine split cutpoint as midpoint of feature fj in this node
  cj <- midpoint(X[[fj]])
  # Partition the data based on the cutpoint
  X_left  <- X[X[[fj]] <= cj, , drop = FALSE]
  X_right <- X[X[[fj]] >  cj, , drop = FALSE]

  # If either side has zero observations, make a leaf (cannot split)
  if (nrow(X_left) == 0 || nrow(X_right) == 0) {
    return(make_leaf(X, D, depth, reason = "empty-child", fj = fj, cj = cj))
  }

  # Evaluate loss for each side with both possible weight assignments (0 or 1)
  loss_left  <- c(loss_fn(0, rownames(X_left),  D), loss_fn(1, rownames(X_left),  D))
  loss_right <- c(loss_fn(0, rownames(X_right), D), loss_fn(1, rownames(X_right), D))

  # Best possible loss on each side
  min_left  <- if (length(loss_left)  > 0) min(loss_left)  else Inf
  min_right <- if (length(loss_right) > 0) min(loss_right) else Inf
  # Combined weighted loss if each side is assigned optimally
  new_loss <- (nrow(X_left) * min_left + nrow(X_right) * min_right) / nrow(X)

  # Decide whether to accept the split
  if (nearly_leq(new_loss, parent_loss)) {


    # Log the accepted split details
    .log("[depth=%d] SPLIT: feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.4f, parent_loss=%.4f",
         depth, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)
    .log("    left losses={%.4f, %.4f}, right losses={%.4f, %.4f}",
         loss_left[1], loss_left[2], loss_right[1], loss_right[2])


    # Determine optimal weight assignment for left and right child nodes
    w_left  <- which.min(loss_left)  - 1  # 0 or 1
    w_right <- which.min(loss_right) - 1
    # Assign weights in D and update local copies X_left, X_right
    D[rownames(X_left),  "w"] <- w_left;  X_left$w  <- w_left
    D[rownames(X_right), "w"] <- w_right; X_right$w <- w_right

    # Recurse deeper, randomizing which side goes first to avoid bias
    if (stats::rbinom(1, 1, 0.5) == 1) {
      left_res  <- split_node(split_feature,        X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature, loss_fn, reduce_weight_fn,
                              max_depth, min_leaf_n, log_fn = log_fn)
      D <- left_res$D  # update global D after left recursion
      right_res <- split_node(split_feature,       X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature, loss_fn, reduce_weight_fn,
                              max_depth, min_leaf_n, log_fn = log_fn)
      D <- right_res$D
    } else {
      right_res <- split_node(split_feature,       X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature, loss_fn, reduce_weight_fn,
                              max_depth, min_leaf_n, log_fn = log_fn)
      D <- right_res$D
      left_res  <- split_node(split_feature,        X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature, loss_fn, reduce_weight_fn,
                              max_depth, min_leaf_n, log_fn = log_fn)
      D <- left_res$D
    }

    # Compute objective at this node after full recursion
    local_obj <- sqrt(sum(D$vsq * D$w, na.rm = TRUE) / (sum(D$w, na.rm = TRUE)^2))
    if (!is.finite(local_obj) || is.nan(local_obj)) {
      local_obj <- Inf
    }

    # Return subtree structure
    return(list(node = fj, split = cj,
                left_tree = left_res, right_tree = right_res,
                `local objective` = local_obj, depth = depth, D = D))
  }

  # If we reach here, the split did not improve the loss (rejected split)
  .log("[depth=%d] REJECTED: feature=%s, cut=%.4f | new_loss=%.16f > parent_loss=%.16f",
       depth, fj, cj, new_loss, parent_loss)

  # Reduce this feature's probability and retry splitting at the same node
  split_feature_updated <- reduce_weight_fn(fj, split_feature)
  return(split_node(split_feature_updated, X, D, parent_loss, depth,
                    explore_proba, choose_feature, loss_fn, reduce_weight_fn,
                    max_depth, min_leaf_n, log_fn = log_fn))
}

#' Fit a single weighted tree for optimized subgroup selection
#'
#' Runs a Double ML procedure to compute pseudo-outcomes, then grows a single decision tree that optimizes an objective function (approximate treatment effect standard error) by reweighting observations. Finally, fits a shallow classifier to summarize the resulting weights.
#'
#' @param data A data frame containing the full dataset (must include outcome, treatment, and sample columns).
#' @param outcome Name of the outcome column.
#' @param treatment Name of the treatment indicator column (0/1).
#' @param sample Name of the sample indicator column (0/1).
#' @param leaf_proba Numeric in `[0,1]`, the probability weight assigned to choosing a "leaf" (no split) at each node during tree building (default 0.25).
#' @param seed Integer seed for reproducibility (default 42). A single `set.seed` is called at the start for all random components.
#' @param verbose Logical, if `TRUE` enables verbose logging of the tree-building process (prints details of splits and decisions). Defaults to `FALSE`.
#' @return A list with components:
#' \item{D}{A data frame (same rows as `data2`) including covariates and additional columns:
#'    `v` (pseudo-outcome), `vsq` (squared pseudo-outcome), `w` (optimized weight, 0/1), and `S` (sample indicator).}
#' \item{f}{An `rpart` model fit on `X` vs `w` (the classifier summarizing the learned weights).}
#' \item{testing_data}{A data frame of the subset of `data` used in tree optimization (specifically, those with `S==1` and finite weights). This is essentially the `data2` from the DML step.}
#' @details **Procedure:** First, cross-fitted pseudo-outcomes are computed using \code{\link{estimate_dml}} (with a default of 5 folds). Then a ridge regression (or GBM, if specified) is used to estimate feature importances for `vsq` (the variance of the pseudo-outcome), which informs the initial feature selection probabilities. A special "leaf" option is included with probability `leaf_proba` to allow the tree to terminate splits. The function then calls an internal recursive splitting routine to assign binary weights `w` to observations, optimizing the treatment effect estimate's precision. The final classifier (`f`) is a decision tree of fixed depth (by default 3) trained on `X` to predict the optimized `w` labels for interpretability.
#'
#' **Randomness and Reproducibility:** The function uses random number generation in the tree-building process (for tie-breaking and exploration). Setting the same `seed` will produce the same result. The `verbose` flag can be used to trace the splitting process for debugging or insight.
#' @examples
#'  sim<-get_data(n=2000,seed=599)
#'  D <-sim$data
#'  dml<-estimate_dml(D,outcome="Yobs",treatment="Tr",sample="S",crossfit= 5)
#'  tr_res <- forest_opt(D,
#'                   outcome="Yobs",treatment="Tr",sample="S",
#'                   leaf_proba=0.25,seed=599)
#' @export
tree_opt <- function(data, outcome, treatment, sample, leaf_proba = 0.25, seed = NULL, verbose = FALSE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(data)) {
      stop(sprintf("Column '%s' not found in data.", col), call. = FALSE)
    }
  }
  check_no_na(data, colnames(data))
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1) {
    stop("`leaf_proba` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("`seed` must be NULL or a single numeric (typically integer) value.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  # Ensure sample indicator has both 0 and 1 present
  if (all(data[[sample]] == 0) || all(data[[sample]] == 1)) {
    stop("Sample indicator column `", sample, "` has no variation (all 0 or all 1). Cannot build tree.", call. = FALSE)
  }


  # 1) Compute cross-fitted pseudo-outcomes
  dml <- estimate_dml(data, outcome = outcome, treatment = treatment, sample = sample, crossfit = 5)
  df_v <- dml$df_v        # data frame of pseudo-outcomes and components
  testing_data <- dml$data2  # corresponding subset of original data (S==1 observations with valid pseudo-outcomes)

  # 2) Extract relevant columns from testing_data
  S_vec  <- testing_data[[sample]]
  Y_vec  <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  X_df   <- testing_data[, setdiff(names(testing_data), c(outcome, treatment, sample)), drop = FALSE]

  # 3) Prepare pseudo-outcome and its variance for feature importance
  v_vals   <- df_v$te
  vsq_vals <- df_v$te_sq

  # 4) Estimate feature importances for vsq:
  if (all(abs(vsq_vals) < .Machine$double.eps) || !is.finite(sum(abs(vsq_vals)))) {
    # Edge case: if vsq is all zero or non-finite, default to uniform feature weights
    feat_prob <- rep(1 / ncol(X_df), ncol(X_df))
    names(feat_prob) <- colnames(X_df)
  } else {
    # Use ridge regression to approximate influence of X on vsq (similar to sklearn's Ridge)
    ridge_fit <- MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1)
    beta <- as.numeric(ridge_fit$coef)  # one coefficient per feature
    names(beta) <- colnames(X_df)
    # Compute feature probabilities proportional to coefficient magnitudes
    if (all(beta == 0) || !is.finite(sum(abs(beta)))) {
      feat_prob <- rep(1 / ncol(X_df), ncol(X_df))
      names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- abs(beta) / sum(abs(beta))
    }
  }
  # Include "leaf" option in feature probabilities
  features <- c("leaf", colnames(X_df))
  proba <- c(leaf_proba, feat_prob)
  proba <- proba / sum(proba)  # normalize just in case
  split_feature <- stats::setNames(proba, features)

  # 5) Build the global state data frame D for tree optimization
  D <- data.frame(X_df, v = v_vals, vsq = vsq_vals, w = rep(1, length(v_vals)), S = S_vec, check.names = FALSE)
  # Use row names as character indices for referencing in `loss` and `split_node`
  rownames(D) <- rownames(X_df) <- as.character(seq_len(nrow(D)))

  # 6) Grow the weighted tree using recursive splitting
  # If verbose, define a logging function to output messages; otherwise, no-op
  log_function <- if (verbose) function(x) cat(x, "\n") else function(x) {}
  run_split <- function() {
    split_node(
      split_feature = split_feature,
      X = D, D = D, parent_loss = Inf, depth = 0L,
      explore_proba = 0.05,
      choose_feature = choose_feature,
      loss_fn = loss,
      reduce_weight_fn = reduce_weight,
      max_depth = 8, min_leaf_n = 5,
      log_fn = log_function
    )
  }
  w_tree <- if (is.null(seed)) run_split() else withr::with_seed(seed, run_split())
  D_final <- w_tree$D  # Updated D with learned weights (`w`)

  # 7) Fit a shallow classifier to summarize the final weights
  final_model <- characterize_tree(X_df, as.integer(D_final$w), max_depth = 3)

  return(list(D = D_final, f = final_model, testing_data = testing_data))
}

#' Ensemble of weighted trees (forest) and Rashomon set selection
#'
#' Builds multiple weighted trees via \code{\link{tree_opt}} logic, then identifies a "Rashomon set" of top-performing trees and aggregates their weight assignments by majority vote.
#'
#' @param data A data frame containing the dataset (must include outcome, treatment, sample indicator).
#' @param outcome Name of the outcome column.
#' @param treatment Name of the treatment indicator column (0/1).
#' @param sample Name of the sample indicator column (0/1).
#' @param leaf_proba Probability mass for the "leaf" option in each tree (see \code{\link{tree_opt}}; default 0.25).
#' @param seed Integer seed for reproducibility (default 42).
#' @param num_trees Number of trees to grow in the forest (default 10).
#' @param vote_threshold Majority vote threshold in (0.5, 1] for assigning final weight=1 (default 2/3, i.e., at least 67% of trees vote 1).
#' @param explore_proba Probability of exploration at leaves in each tree (default 0.05, see \code{\link{split_node}}).
#' @param feature_est Method for feature importance estimation, `"Ridge"` (default) or `"GBM"`.
#' @param top_k_trees Logical. If `TRUE`, selects the Rashomon set as the top `k` trees with lowest objective. If `FALSE`, uses a cutoff threshold.
#' @param k Integer, number of top trees to keep if `top_k_trees = TRUE` (default 10).
#' @param cutoff Either `"baseline"` or a numeric value. If `top_k_trees = FALSE`, this is the loss cutoff for selecting Rashomon set (if `"baseline"`, uses the baseline loss of not weighting any observations).
#' @param verbose Logical, if `TRUE`, prints progress messages (e.g., estimated PATE and number of trees selected). Default `FALSE`.
#' @return A list with components:
#' \item{D_rash}{A data frame similar to `D_forest` but containing only the Rashomon set of trees' weight columns and additional columns: `w_opt` (the ensemble-voted weight for each observation) and `vote_count` (how many trees voted to include that observation).}
#' \item{D_forest}{A data frame with all trees' weight assignments. It contains covariates, `v`, `vsq`, `S`, `lX` (overlap factor = 1/b), and columns `w_tree_1, ..., w_tree_num_trees` for each tree's weights.}
#' \item{w_forest}{A list of length `num_trees`, where each element is the full tree object returned by the internal split routine (containing the tree structure and local objective).}
#' \item{rashomon_set}{An integer vector of indices of trees that were selected into the Rashomon set.}
#' \item{f}{An `rpart` model fit on `X` vs `w_opt` (the final classifier summarizing the ensemble's selected subgroup).}
#' \item{testing_data}{The data frame of observations used in tree construction (same as in `tree_opt` and output from `estimate_dml`, i.e., those with `S==1` and finite pseudo-outcomes).}
#' @details Each tree in the forest is built on the same cross-fitted pseudo-outcomes (from one call to \code{\link{estimate_dml}}) for consistency. If `feature_est = "GBM"`, a gradient boosting model (from \pkg{gbm}) is used instead of ridge regression to estimate feature importance for splitting probabilities. The Rashomon set is defined as either the top `k` trees with lowest objective (if `top_k_trees=TRUE`), or all trees with objective below a cutoff. The baseline loss is defined as the loss with no selection (all weights = 1 for S==1 group).
#' @examples
#'  sim<-get_data(n=2000,seed=599)
#'  D <-sim$data
#'  dml<-estimate_dml(D,outcome="Yobs",treatment="Tr",sample="S",crossfit= 5)
#'  fo_res<-forest_opt(D,
#'                   outcome="Yobs",treatment="Tr",sample="S",
#'                   leaf_proba=0.25,seed=3,
#'                   num_trees=50, vote_threshold=2/3,
#'                   explore_proba=0.05,feature_est="Ridge",
#'                   top_k_trees=FALSE,cutoff="baseline")
#' @export
forest_opt <- function(data,
                       outcome,
                       treatment,
                       sample,
                       leaf_proba = 0.25,
                       seed = NULL,
                       num_trees = 10,
                       vote_threshold = 2/3,
                       explore_proba = 0.05,
                       feature_est = "Ridge",   # "Ridge" or "GBM"
                       top_k_trees = FALSE,
                       k = 10,
                       cutoff = "baseline",
                       verbose = FALSE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(data)) {
      stop(sprintf("Column '%s' not found in data.", col), call. = FALSE)
    }
  }
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1) {
    stop("`leaf_proba` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(num_trees) || num_trees < 1) {
    stop("`num_trees` must be a positive integer.", call. = FALSE)
  }
  num_trees <- as.integer(num_trees)
  if (!is.numeric(vote_threshold) || vote_threshold <= 0 || vote_threshold > 1) {
    stop("`vote_threshold` must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(explore_proba) || explore_proba < 0 || explore_proba > 1) {
    stop("`explore_proba` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.character(feature_est) || length(feature_est) != 1 || !(tolower(feature_est) %in% c("ridge", "gbm"))) {
    stop("`feature_est` must be \"Ridge\" or \"GBM\" (string).", call. = FALSE)
  }
  if (!is.logical(top_k_trees) || length(top_k_trees) != 1) {
    stop("`top_k_trees` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(k) || k < 1) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  if (! (is.numeric(cutoff) || (is.character(cutoff) && cutoff == "baseline")) ) {
    stop("`cutoff` must be \"baseline\" or a numeric value.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (all(data[[sample]] == 0) || all(data[[sample]] == 1)) {
    stop("Sample indicator column `", sample, "` has no variation (all 0 or all 1). Cannot build forest.", call. = FALSE)
  }

  # Set random seed for reproducibility
  set.seed(seed)

  # 1) Compute cross-fitted pseudo-outcomes (fixed 5-fold by default for forest)
  dml <- estimate_dml(data, outcome = outcome, treatment = treatment, sample = sample, crossfit = 5)
  df_v <- dml$df_v
  testing_data <- dml$data2

  # 2) Prepare data and columns
  S_vec <- testing_data[[sample]]
  Y_vec <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  X_df <- testing_data[, setdiff(names(testing_data), c(outcome, treatment, sample)), drop = FALSE]
  n  <- nrow(testing_data)
  v_vals   <- df_v$te
  vsq_vals <- df_v$te_sq

  if (verbose) {
    message(sprintf("PATE Estimate (mean of v): %.4f", mean(v_vals)))
  }

  # 3) Determine feature selection probabilities for splitting
  features <- c("leaf", colnames(X_df))
  if (tolower(feature_est) == "ridge") {
    if (all(abs(vsq_vals) < .Machine$double.eps) || !is.finite(sum(vsq_vals))) {
      feat_prob <- rep(1 / ncol(X_df), ncol(X_df))
      names(feat_prob) <- colnames(X_df)
    } else {
      ridge_fit <- MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1)
      beta <- as.numeric(ridge_fit$coef); names(beta) <- colnames(X_df)
      if (all(beta == 0) || !is.finite(sum(abs(beta)))) {
        feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
      } else {
        feat_prob <- abs(beta) / sum(abs(beta))
      }
    }
  } else {  # "GBM"
    df_gbm <- data.frame(X_df, vsq = vsq_vals)
    fit_gbm <- function() {
      gbm::gbm(
        vsq ~ ., data = df_gbm, distribution = "gaussian",
        n.trees = 100, interaction.depth = 3, shrinkage = 0.1,
        bag.fraction = 0.5, train.fraction = 1.0, n.minobsinnode = 2, verbose = FALSE
      )
    }
    gbm_fit <- if (is.null(seed)) fit_gbm() else withr::with_seed(seed + 1000L, fit_gbm())
    rel_inf <- gbm::relative.influence(gbm_fit, n.trees = 100, sort. = FALSE)
    names(rel_inf) <- colnames(X_df)
    if (sum(rel_inf) <= 0 || !is.finite(sum(rel_inf))) {
      feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- rel_inf / sum(rel_inf)
    }
  }
  proba <- c(leaf_proba, feat_prob)
  proba <- proba / sum(proba)
  split_feature <- stats::setNames(proba, features)

  # 4) Build forest of weighted trees
  # Initialize a data frame to collect weights from all trees
  D_forest <- data.frame(X_df, v = v_vals, vsq = vsq_vals, S = S_vec, stringsAsFactors = FALSE)
  # Overlap factor lX = 1/b for each observation (for diagnostics, from df_v$b)
  D_forest$lX <- 1 / df_v$b
  # Storage for individual tree results
  w_forest <- vector("list", num_trees)

  for (t_idx in seq_len(num_trees)) {
    # Initialize all weights = 1 for this tree
    D_tree <- data.frame(X_df, v = v_vals, vsq = vsq_vals, w = rep(1, n), S = S_vec, stringsAsFactors = FALSE)
    rownames(D_tree) <- as.character(seq_len(n))
    rownames(X_df) <- as.character(seq_len(n))
    # Grow one weighted tree
    grow_one <- function() {
      split_node(
        split_feature = split_feature,
        X = D_tree, D = D_tree, parent_loss = Inf, depth = 0L,
        explore_proba = explore_proba,
        choose_feature = choose_feature,
        loss_fn = loss,
        reduce_weight_fn = reduce_weight
      )
    }
    tree_res <- if (is.null(seed)) grow_one() else withr::with_seed(seed + t_idx, grow_one())
    D_updated <- tree_res$D
    # Store weights from this tree and the tree object
    D_forest[[paste0("w_tree_", t_idx)]] <- D_updated$w
    w_forest[[t_idx]] <- tree_res
  }

  # 5) Select Rashomon set of trees based on objective values
  obj_values <- vapply(w_forest, function(res) res[["local objective"]], numeric(1))
  tree_indices <- seq_along(w_forest)
  if (top_k_trees) {
    if (k > num_trees) {
      warning("k (", k, ") is greater than num_trees; using k = num_trees.", call. = FALSE)
      k <- num_trees
    }
    ord <- order(obj_values, na.last = NA)  # ascending order (lower is better)
    rashomon_set <- utils::head(ord, k)
  } else {
    if (identical(cutoff, "baseline")) {
      # baseline loss: no selection (all w for S==1 group = 1)
      baseline_loss <- sqrt(sum(D_forest$vsq) / (nrow(D_forest)^2))
      cutoff_val <- baseline_loss
    } else {
      cutoff_val <- as.numeric(cutoff)
    }
    rashomon_set <- tree_indices[ obj_values < cutoff_val ]
  }
  if (length(rashomon_set) == 0) {
    warning("No trees selected into Rashomon set (all trees have objective above cutoff).", call. = FALSE)
  }
  not_in_set <- setdiff(tree_indices, rashomon_set)

  # 6) Construct datasets for selected trees
  D_rash <- D_forest
  if (length(not_in_set) > 0) {
    drop_cols <- paste0("w_tree_", not_in_set)
    keep_cols <- setdiff(names(D_rash), drop_cols)
    D_rash <- D_rash[, keep_cols, drop = FALSE]
  }
  # Subset of weight columns for rashomon set (if any)
  weight_cols <- grep("^w_tree_", names(D_rash), value = TRUE)
  D_weights <- if (length(weight_cols) > 0) D_rash[, weight_cols, drop = FALSE] else data.frame()

  # 7) Majority vote ensemble of selected trees
  if (ncol(D_weights) > 0) {
    row_means <- rowMeans(D_weights)
    D_rash$w_opt <- as.integer(row_means > vote_threshold)
    D_rash$vote_count <- rowSums(D_weights)
  } else {
    D_rash$w_opt <- integer(nrow(D_rash))  # all zeros if no tree in set
    D_rash$vote_count <- integer(nrow(D_rash))
  }

  # 8) Final classifier on ensemble weights
  final_classifier <- characterize_tree(X_df, as.factor(D_rash$w_opt))

  return(list(D_rash = D_rash,
              D_forest = D_forest,
              w_forest = w_forest,
              rashomon_set = rashomon_set,
              f = final_classifier,
              testing_data = testing_data))
}




