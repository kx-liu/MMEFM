# MMEFM

MMEFM fits the Multilevel Main Effects Matrix Factor Model to groups of matrix-valued time series. Groups share observation times but may have different row and column dimensions. The model separates global and local grand means, row main effects, column main effects, and bilinear common components.

This separation distinguishes changes in overall levels, rows, or columns from row-column interactions. Global structures describe shared dynamics; local structures retain group-specific variation. Main effects may be nonstationary and remain estimable at individual time points under the manuscript's theoretical conditions.

## Installation

MMEFM is under development and is not yet on CRAN. Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("kx-liu/MMEFM")
library(MMEFM)
```

## Quick start

Simulate three named groups with unequal matrix dimensions, then fit the model at the known generating ranks. This is the lightweight example used in the vignette.

```r
library(MMEFM)

rank <- list(r1 = 1L, l1 = 1L, r2 = c(a = 1L, b = 1L, c = 1L), l2 = c(a = 1L, b = 1L, c = 1L),
             kr = 1L, kc = 1L, kr_m = c(a = 1L, b = 1L, c = 1L), kc_m = c(a = 1L, b = 1L, c = 1L))
simulation <- gen_MMEFM(T = 40L, p = c(a = 6L, b = 7L, c = 8L), q = c(a = 7L, b = 8L, c = 6L),
                        rank = rank, burn = 20L, filter_length = 40L, seed = 2026L)
fit <- est_MMEFM(simulation$Xt, rank = simulation$rank, K0 = 8L, max_iter = 50L, tol = 1e-4, seed = 2026L)
summary(fit)
predicted <- fitted(fit)
global_common <- fitted(fit, component = "global_common")
errors <- residuals(fit)
lapply(predicted, dim)
```

The summary reports dimensions, ranks, and method choices. To illustrate global common-factor detection and group screening, use a small serial bootstrap:

```r
detection <- detect_MMEFM_global(fit, B = 3L, K0 = 8L, max_iter = 50L, tol = 1e-4,
                                 parallel = FALSE, seed = 2026L)
detection
detection$selected
```

Three bootstrap replications do not provide reliable statistical calibration. The example is not guaranteed to detect all active groups.

## Main functionality

- **Simulation:** `gen_MMEFM()` generates grouped data and model truth. `gen_piecewise_arfima()` generates scalar series with a stable AR component and segment-specific fractional memory. Both support Gaussian or raw Student-t innovations.
- **Rank selection:** `select_MMEFM_rank()` estimates eight global/local loading rank classes using unperturbed consecutive eigenvalue ratios. It checks the feasibility of the selected ranks and errors on infeasible combinations without repair. A small noisy sample need not yield the true ranks.
- **Estimation:** `est_MMEFM()` accepts known ranks or selects them automatically. The resulting `mmefm_fit` supports `summary()`, `fitted()`, and `residuals()`.
- **Detection:** `detect_MMEFM_global()` assesses whether a global common component exists and screens the groups sharing it. The second-largest group statistic and one circular-shift bootstrap cutoff are used for existence detection and screening.

## Interpretation

The input `Xt` is a list of numeric arrays with dimensions `T x p_m x q_m`, ordered by time, rows, and columns. The final fitted signal combines estimated low-rank global/local main effects, including refitted global loadings, with global/local common components. Direct moment estimates initialize the procedure; they are not substituted for the final main-effect factor representations.

`fit$check_Y` is the intermediate residual after removing direct additive estimates, not the final residual returned by `residuals(fit)`. Convergence diagnostics are in `fit$convergence`, not `summary(fit)`.

Whole-sample circular shifts require a single stationary segment. They are not automatically justified for piecewise nonstationary data. Assess ranks and convergence before interpreting fitted components; a successful fit alone does not establish the model assumptions.

## Further documentation

The [Getting Started vignette](vignettes/getting-started.Rmd) explains the model, estimation stages, rank conventions, and worked example. Use `?gen_piecewise_arfima`, `?gen_MMEFM`, `?select_MMEFM_rank`, `?est_MMEFM`, and `?detect_MMEFM_global` for complete function documentation.

## Accompanying manuscript

Kaixin Liu, Clifford Lam, and Zetai Cen. *Multilevel Main Effects Matrix Factor Model*. Manuscript.

## Development status and license

This is an unreleased development package. MMEFM is licensed under the MIT license; see [LICENSE](LICENSE).
