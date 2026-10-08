# MMEFM

MMEFM fits the Multilevel Main Effects Matrix Factor Model to multiple matrix-valued time series. It preserves matrix structure while separating shared and group-specific effects: global and local grand means, row main effects, column main effects, and matrix-factor interactions. The model allows nonstationary main effects whose time-specific values can still be consistently estimated under the manuscript's theoretical conditions.

## Installation

MMEFM is under development and is not yet on CRAN. Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("kx-liu/MMEFM")
library(MMEFM)
```

## Quick start

This example uses the known simulation ranks. The detection call uses **B = 3 only for demonstration**; three draws do not provide reliable production bootstrap calibration.

```r
library(MMEFM)

rank <- list(r1 = 1L, l1 = 1L, r2 = c(a = 1L, b = 1L, c = 1L), l2 = c(a = 1L, b = 1L, c = 1L),
             kr = 1L, kc = 1L, kr_m = c(a = 1L, b = 1L, c = 1L), kc_m = c(a = 1L, b = 1L, c = 1L))
simulation <- gen_MMEFM(T = 40L, p = c(a = 6L, b = 7L, c = 8L), q = c(a = 7L, b = 8L, c = 6L),
                        rank = rank, burn = 20L, filter_length = 40L, seed = 2026L)
fit <- est_MMEFM(simulation$Xt, rank = simulation$rank, K0 = 8L, max_iter = 50L, seed = 2026L)
summary(fit)
predicted <- fitted(fit)
errors <- residuals(fit)
lapply(predicted, dim)

detection <- detect_MMEFM_global(fit, B = 3L, K0 = 8L, max_iter = 50L, parallel = FALSE, seed = 2026L)
detection
```

See the [Getting Started vignette](vignettes/getting-started.Rmd) for rank conventions, component reconstruction, and interpretation.

## Main functionality

- `gen_MMEFM()` simulates heterogeneous groups with Gaussian or raw Student-t innovations, general piecewise-memory specifications, optional global ownership sets, and optional full component storage.
- `select_MMEFM_rank()` applies the manuscript's eigenvalue-ratio rules to select eight rank classes independently, then checks combined-rank feasibility. Infeasible combinations error rather than being repaired.
- `est_MMEFM()` supports supplied and automatic ranks and returns an `mmefm_fit`. Use `summary()`, `fitted()` (including individual components), and `residuals()` to inspect the fit.
- `detect_MMEFM_global()` uses one bootstrap cutoff for existence detection and group screening. Serial and PSOCK parallel bootstrap are available; process startup may make parallel execution slower for small B. When running Monte Carlo replications in parallel externally, normally leave detector-level parallelization off.

Use `?gen_MMEFM`, `?select_MMEFM_rank`, `?est_MMEFM`, and `?detect_MMEFM_global` for argument and return-value details.

## Statistical conventions

`Xt` is a list of finite numeric arrays with dimensions `T x p_m x q_m`: groups share T but may have different row and column dimensions. Model ranks are currently positive. Under IC1 centering, each global/local row rank sum is at most `p_m - 1`, and each column sum is at most `q_m - 1`; equality is allowed. Invalid ranks and singular unregularized numerical systems raise errors without silent repair.

Simulation ownership (`active_global`) determines which groups carry global common interactions; it does not change nominal positive global ranks. Whole-sample circular-shift detection assumes a stationary segment and is not justified for arbitrary nonstationary piecewise series merely because the generator supports them.

## Accompanying manuscript

The methodology implemented in MMEFM is developed in:

Kaixin Liu, Clifford Lam, and Zetai Cen. *Multilevel Main Effects Matrix Factor Model*. Manuscript.

## Development status and license

This is an unreleased development package. MMEFM is licensed under the MIT license; see [LICENSE](LICENSE).
