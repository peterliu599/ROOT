## Package Updates

### Changes in Version 0.2.1

#### Bug fixes (generalizability mode)

* `compute_transport_scores()` now includes the overall selection-odds
  prefactor pi/(1 - pi) = n1/n0 in the IPW score `v`, as in Parikh et al.
  (2025). Previously the reported unweighted and weighted estimates and their
  SEs were inflated by (1 - pi)/pi.
* New default objective `objective_generalizability()`: the finite-sample SE
  of the IPW WTATE estimator, with residuals centered at the mean of the kept
  units and scaled by `kappa_w`, which re-normalizes the selection odds to
  pi_w = P(S = 1 | w(X) = 1). Previously the search minimized an uncentered
  second moment, sqrt(sum(w * v^2)) / sum(w), which can select different trees.
* The weighted estimate and SE are now `kappa_w * mean(v[A])` and
  `kappa_w * sd(v[A]) / sqrt(n_A)`; `kappa_w` is returned in `$estimate`.

#### Terminology

* The unweighted generalizability estimate is now labelled `"STATE"` (sample
  target average treatment effect, following Parikh et al. 2025) instead of
  `"SATE"`, which conventionally denotes the trial-sample average treatment
  effect. The weighted estimate remains `"WTATE"`.
* TATE and WTATE are now described as *target* (not *transported*) average
  treatment effects, and `v` is described as an inverse propensity weighting
  (IPW) transported score.
* `summary()` and `print()` now label the estimates block
  "Estimand summary (generalizability mode)".

#### Bug fixes (general mode)

* When only a `v` column is supplied, `vsq` is no longer precomputed around
  the full-sample mean; `objective_default()` now centers `v` at the mean of
  the currently included units at every evaluation.
* When the data contain neither `vsq` nor `v` and no `global_objective_fn`
  is supplied, ROOT uses `vsq = 1` and now issues a message explaining that
  the default objective then favors keeping all observations.

#### Documentation

* The portfolio example (README, `optimization_path_example` vignette, and
  `?ROOT`) now minimizes the variance of the equally weighted portfolio,
  w' Sigma w / (1'w)^2, using the full return covariance matrix, instead of
  summing marginal variances.
* The `generalizability_path_example` vignette interpretation now reflects
  the corrected results and checks them against the data-generating process.
* The quickstart vignette now correctly describes how the default objective
  reads the `vsq` or `v` column (it is not the first column of `data`).

### Changes in Version 0.2.0 (2026-04-09)

#### New features
*  Added `positivity_trim` argument to `ROOT()` and `characterizing_underrep()` for pre-trimming target observations with structural positivity violations prior to optimization. Supports a fixed numeric threshold (e.g., `positivity_trim = 0.05`) or the data-driven rule of Crump, Hotz, Imbens & Mitnik (2009, *Biometrika*) via `positivity_trim = "crump"`. Only active in generalizability mode. Trimming metadata is stored in the returned object as `positivity_trim_info` and reported by `summary()`.
*  Added `split_strategy` argument to `ROOT()` and `characterizing_underrep()` with options `"midpoint"` (default, backward-compatible) and `"random_quantile"`. The random quantile strategy draws a uniform random quantile within `quantile_range` (default 10th–90th percentile) at each split, improving covariate-space coverage across the forest. The strategy used is stored in the returned object and reported by `summary()`.
*  New internal helpers: `crump_alpha()`, `trim_positivity_violations()`, and `split_point()`.

#### Documentation
*  Added roxygen documentation for all new arguments and helper functions.
*  Added Crump et al. (2009) to `@references` in `ROOT()`.
*  `summary.ROOT()` now reports the split strategy and positivity trimming diagnostics.

#### Testing
*  New test files: `test-positivity.R` and `test-split-strategy.R`.

### Changes in Version 0.1.1 (2026-03-04)
*  Updated help files
*  Expanded unit testing
*  Updated authorship

### Changes in Version 0.1.0 (2026-02-23)
*  First version released on GitHub (<https://github.com/peterliu599/ROOT>)
