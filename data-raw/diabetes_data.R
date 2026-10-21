# =============================================================================
# diabetes_data: simulated trial + target data for the generalizability example
# =============================================================================
#
# Setting (illustrative). A randomized trial of a new glucose-lowering drug,
# with outcome = reduction in HbA1c (percentage points) at 6 months. We want
# to generalize the trial result to a broader target population of adults with
# type 2 diabetes.
#
# The design deliberately contains all four combinations of
# (under-sampled in the trial?) x (modifies the treatment effect?):
#
#                    modifies effect                    does not modify effect
#   under-sampled    CKD, HeartFailure, Insulin,        Race_Black, Rural
#                    Age >= 70
#   well sampled     HbA1c (mild, linear)               Sex_Male, BMI, DietYes
#
# so that the example shows ROOT flagging the under-sampled effect modifiers
# while keeping under-sampled groups whose treatment effect is no different.
#
# Data-generating process (one pool of n = 10,000 people):
#
#   Covariates
#     Age          ~ N(58, 12^2), truncated to [30, 85], rounded to whole years
#     Sex_Male     ~ Bernoulli(0.50)
#     Race_Black   ~ Bernoulli(0.25)
#     Rural        ~ Bernoulli(0.20)
#     BMI          ~ N(31, 5^2), truncated to [18, 50], 1 decimal
#     HbA1c        ~ N(8.2, 1.1^2), truncated to [6.5, 12], 1 decimal (baseline, %)
#     Insulin      ~ Bernoulli(expit(-1.4 + 0.04 (Age - 58) + 0.5 (HbA1c - 8)))
#     CKD          ~ Bernoulli(expit(-2.0 + 0.06 (Age - 58)))   # chronic kidney disease
#     HeartFailure ~ Bernoulli(expit(-2.0 + 0.03 (Age - 58)))
#     DietYes      ~ Bernoulli(0.15)
#
#   Trial participation (S = 1 trial, S = 0 target)
#     logit P(S = 1 | X) = -0.9 - 1.25 CKD - 0.8 HeartFailure - 0.4 Insulin
#                          - 0.018 (Age - 58) - 0.5 Race_Black - 0.35 Rural
#     (about 2,000 trial participants and 8,000 target units)
#
#   Treatment (trial only): Tr ~ Bernoulli(0.5)
#
#   Conditional average treatment effect
#     tau(X) = 1 + 0.1 (HbA1c - 8) - 2.0 CKD - 2.4 HeartFailure + 1.5 Insulin
#              + 2.0 1{Age >= 70}
#
#   Outcome (trial only)
#     Y = mu0(X) + Tr * tau(X) + eps,   eps ~ N(0, 0.5^2)
#     mu0(X) = 0.2 (HbA1c - 8) + 0.1 Insulin - 0.2 Race_Black + 0.1 Rural
#              - 0.01 (BMI - 31) + 0.2 DietYes + 0.1 HeartFailure
#
# Notes
#   * The effects are illustrative: the drug is counterproductive for people
#     with CKD or heart failure and much more effective for insulin users and
#     people aged 70+.
#   * Race_Black and Rural shift the outcome level (prognostic) but not the
#     treatment effect.
#   * The selection model is linear-logistic in the covariates, so the
#     main-effects logistic regression used by ROOT is correctly specified.
#   * True CATEs are NOT stored in the data (ROOT would otherwise use them as
#     splitting covariates); use true_cate() below to compute them.
# =============================================================================

simulate_diabetes_data <- function(n = 10000, seed = 2026) {
  set.seed(seed)

  # --- 1) Covariates --------------------------------------------------------
  Age          <- round(pmin(pmax(rnorm(n, 58, 12), 30), 85))
  Sex_Male     <- rbinom(n, 1, 0.50)
  Race_Black   <- rbinom(n, 1, 0.25)
  Rural        <- rbinom(n, 1, 0.20)
  BMI          <- round(pmin(pmax(rnorm(n, 31, 5), 18), 50), 1)
  HbA1c        <- round(pmin(pmax(rnorm(n, 8.2, 1.1), 6.5), 12), 1)
  Insulin      <- rbinom(n, 1, plogis(-1.4 + 0.04 * (Age - 58) + 0.5 * (HbA1c - 8)))
  CKD          <- rbinom(n, 1, plogis(-2.0 + 0.06 * (Age - 58)))
  HeartFailure <- rbinom(n, 1, plogis(-2.0 + 0.03 * (Age - 58)))
  DietYes      <- rbinom(n, 1, 0.15)

  # --- 2) Trial participation ----------------------------------------------
  pi_true <- plogis(-0.9 - 1.25 * CKD - 0.8 * HeartFailure - 0.4 * Insulin -
                      0.018 * (Age - 58) - 0.5 * Race_Black - 0.35 * Rural)
  S <- rbinom(n, 1, pi_true)

  # --- 3) Treatment, randomized within the trial only ------------------------
  Tr <- rep(NA_integer_, n)
  Tr[S == 1] <- rbinom(sum(S == 1), 1, 0.5)

  # --- 4) Outcomes -----------------------------------------------------------
  tau <- true_cate(Age = Age, HbA1c = HbA1c, CKD = CKD,
                   HeartFailure = HeartFailure, Insulin = Insulin)
  mu0 <- 0.2 * (HbA1c - 8) + 0.1 * Insulin - 0.2 * Race_Black +
    0.1 * Rural - 0.01 * (BMI - 31) + 0.2 * DietYes + 0.1 * HeartFailure
  eps <- rnorm(n, 0, 0.5)
  Y <- rep(NA_real_, n)
  Y[S == 1] <- mu0[S == 1] + Tr[S == 1] * tau[S == 1] + eps[S == 1]

  data.frame(
    Age, Sex_Male, Race_Black, Rural, BMI, HbA1c, Insulin, CKD, HeartFailure,
    DietYes, S = as.integer(S), Tr = as.integer(Tr), Y
  )
}

# True conditional average treatment effect (for checking results)
true_cate <- function(Age, HbA1c, CKD, HeartFailure, Insulin) {
  1 + 0.1 * (HbA1c - 8) - 2.0 * CKD - 2.4 * HeartFailure + 1.5 * Insulin +
    2.0 * (Age >= 70)
}

# ---- Generate, check, and save ----------------------------------------------
diabetes_data <- simulate_diabetes_data(n = 10000, seed = 2026)

# Quick checks: sample sizes, covariate balance, and true estimands
table(S = diabetes_data$S)
round(t(sapply(diabetes_data[, 1:10], function(x) tapply(x, diabetes_data$S, mean))), 3)
with(subset(diabetes_data, S == 0), {
  tau  <- true_cate(Age, HbA1c, CKD, HeartFailure, Insulin)
  keep <- CKD == 0 & HeartFailure == 0 & Insulin == 0 & Age < 70
  c(TATE_true = mean(tau), WTATE_true = mean(tau[keep]), share_kept = mean(keep))
})

usethis::use_data(diabetes_data, overwrite = TRUE)
