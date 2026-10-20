ROOT: Weighted Trees/Forests for ROOT-style Functional Optimization
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/peterliu599/ROOT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/peterliu599/ROOT/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/peterliu599/ROOT/graph/badge.svg)](https://app.codecov.io/gh/peterliu599/ROOT)
<!-- badges: end -->

> Rashomon-set Optimal Trees for interpretable functional optimization
> and treatment effect generalizability

------------------------------------------------------------------------

# Overview

**ROOT (Rashomon set of Optimal Trees)** is a **general functional
optimization framework** for learning interpretable binary weight
functions, represented as sparse decision trees.

At its core, ROOT searches over functions $w(X) \in \{0,1\}$ that assign
inclusion or exclusion weights to units based on their covariates $X$.
The optimization objective (or **loss function**) can be chosen to
reflect different scientific goals.

- **Framework flexibility:** By substituting different loss functions,
  ROOT can be applied to diverse settings.
- **Interpretability:** Solutions are restricted to sparse decision
  trees, ensuring outputs are human-readable and suitable for
  communication with collaborators.
- **Rashomon principle:** Rather than return a single solution, ROOT
  identifies a **Rashomon set** of near-optimal trees, highlighting
  multiple plausible explanations. A **characteristic tree** can then be
  extracted to summarize common patterns.

------------------------------------------------------------------------

## ROOT for generalizability

The original ROOT paper demonstrates this framework in the context of
**generalizing causal effects from randomized trials to a target
population**.

- **Problem:** Randomized trials are internally valid but often differ
  from the target population. If trial participants underrepresent
  certain effect-modifying subgroups, estimates of the **Target Average
  Treatment Effect (TATE)** can be imprecise or misleading.
- **Solution:** ROOT introduces a refined estimand, the **Weighted
  Target Average Treatment Effect (WTATE)**, defined over the subset of
  the population that is sufficiently represented.
  - Units are either included ($w(X) = 1$) or excluded ($w(X) = 0$)
    according to learned tree-structured rules.
  - The loss function is designed to minimize the variance of the
    estimator while retaining as much of the target population as
    possible.
- **Interpretation:** ROOT’s decision trees explicitly describe *who is
  underrepresented* and *why they are excluded*. The Rashomon set
  provides multiple, equally valid ways of describing these subgroups.

For a detailed worked example of ROOT in generalizability mode, see the
[`generalizability_path_example`
vignette](vignettes/generalizability_path_example.Rmd).

------------------------------------------------------------------------

## Installation

You can install the development version of `ROOT` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("peterliu599/ROOT")
```

------------------------------------------------------------------------

## Example: Portfolio selection via variance minimization

ROOT can be applied to any optimization problem that can be expressed as
a binary inclusion/exclusion decision. In this example, we use ROOT to
choose which of 100 simulated assets to include in an **equally weighted,
low-variance portfolio**. Each asset is characterized by its market beta
and annualized volatility, and ROOT learns an interpretable tree-structured
rule describing which assets to include ($w = 1$).

The objective is the return variance of the equally weighted portfolio of
the included assets, $L(w) = w^\top \Sigma w / (\mathbf{1}^\top w)^2$, where
$\Sigma$ is the return covariance matrix. Because every pair of included
assets contributes a covariance term, the objective cannot be split into
independent per-asset contributions; it is supplied to ROOT as a custom
`global_objective_fn`.

``` r
library(ROOT)
set.seed(123)

n_assets <- 100

# Asset features
volatility <- runif(n_assets, 0.05, 0.40)  # annualised volatility
beta       <- runif(n_assets, 0.5,  1.8)   # market beta
sector     <- sample(c("Tech", "Finance", "Energy", "Health"),
                     n_assets, replace = TRUE)

# Simulate returns from a single-factor market model
market      <- rnorm(1000, 0.0005, 0.01)
returns_mat <- sapply(seq_len(n_assets), function(i)
  beta[i] * market + rnorm(1000, 0, volatility[i] / sqrt(252))
)
Sigma_hat <- cov(returns_mat)  # estimated return covariance matrix

# Variance of the equally weighted portfolio of included assets
portfolio_variance <- function(D) {
  w <- D$w
  if (sum(w) == 0) return(Inf)
  sum(Sigma_hat[w == 1, w == 1]) / sum(w)^2  # = w' Sigma w / (1'w)^2
}

dat_portfolio <- data.frame(
  vsq    = diag(Sigma_hat),  # marginal variances; guide split sampling only
  vol    = volatility,
  beta   = beta,
  sector = as.integer(factor(sector))
)

portfolio_fit <- ROOT(
  data                = dat_portfolio,
  global_objective_fn = portfolio_variance,
  num_trees           = 20,
  top_k_trees         = TRUE,
  k                   = 10,
  seed                = 42
)

summary(portfolio_fit)  # full summary
print(portfolio_fit)    # brief summary
plot(portfolio_fit)     # characterized tree
```

The characterized tree gives a plain-language rule for which asset types
belong in the portfolio, and `portfolio_fit$D_rash$w_opt` holds the final
inclusion decision for each asset. In this example the rule is driven mainly
by market beta: ROOT keeps 33 mostly low-beta assets (beta below about 0.98),
and the equally weighted portfolio of those assets has about 47% lower return
variance than the equally weighted portfolio of all 100 assets.

For a detailed worked example of ROOT in optimization mode, see the
[`optimization_path_example`
vignette](vignettes/optimization_path_example.Rmd).

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/peterliu599/ROOT/issues).

------------------------------------------------------------------------

## License

The contents of this repository are distributed under the MIT license.

------------------------------------------------------------------------

## Reference

Parikh, H., Ross, R. K., Stuart, E., & Rudolph, K. E. (2025). Who Are We
Missing?: A Principled Approach to Characterizing the Underrepresented
Population. *Journal of the American Statistical Association*, 120(551),
1414–1423. <https://doi.org/10.1080/01621459.2025.2495319>
