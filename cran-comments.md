## Resubmission

This is a resubmission following CRAN's manual review of MMEFM 0.1.0. MMEFM 0.1.1 is intended to be the first CRAN release; no version has yet been accepted on CRAN.

In response to the reviewer's request, the time-length argument `T` has been renamed to `TT` in `gen_MMEFM()` and `gen_piecewise_arfima()`. An internal standalone variable `F` has also been renamed to `F_local`. All affected documentation, examples, tests, and vignettes have been updated.

The statistical methods, numerical behavior, and returned object schemas remain unchanged. Deterministic before/after comparisons produced identical full generator outputs and ARFIMA series, with caller RNG state and kind preserved.

## Test environments

* Local: macOS 27.0.1, R 4.6.1 (2026-06-24), aarch64-apple-darwin23; full `R CMD check --as-cran`.
* GitHub Actions uses `R CMD check --no-manual` on macOS and Windows with R release, and Ubuntu with R release and R-devel. These checks are distinct from the full local check reported below.

## R CMD check results

`R CMD check --as-cran` on `MMEFM_0.1.1.tar.gz`: 0 ERRORs, 0 WARNINGs, 1 NOTE.

* CRAN incoming feasibility: "New submission". This is expected because MMEFM has not previously been released on CRAN.

HTML validation, including mathematical rendering, completed successfully with V8 8.2.0 and HTML Tidy 5.8.0. Evaluated examples, vignette checks and rebuilding, and the PDF manual passed. All 2,324 regression-test assertions passed, with 0 failures, warnings, or skips.

## Additional remarks

MMEFM implements statistical methodology for grouped matrix-valued time series, separating global and local main effects from bilinear common components. The word "nonstationary" is an established technical term in time series analysis.
