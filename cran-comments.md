## Release status

First submission of MMEFM, version 0.1.0. This is a release candidate; the package has not been submitted to CRAN.

## Test environments

* Local: R 4.6.1 (2026-06-24), aarch64-apple-darwin23, macOS 27.0.1.
* Final release-candidate Windows, Linux, and R-devel checks are pending.

## R CMD check results

`R CMD check --as-cran` on `MMEFM_0.1.0.tar.gz`: 0 ERRORs, 0 WARNINGs, 2 NOTEs.

* CRAN incoming feasibility: "New submission". This is expected for the first submission; no development-version NOTE remains.
* HTML manual: "Skipping checking math rendering: package 'V8' unavailable". HTML Tidy 5.8.0 was selected through `R_TIDYCMD` and HTML validation reported no problems, but mathematical HTML rendering was not validated. Full mathematical HTML validation is pending.

All 2,191 regression-test assertions, evaluated examples, vignette checks, and the PDF manual passed. The source build took 3.86 seconds and the full local check took 29.66 seconds.

## Additional remarks

MMEFM implements statistical methodology for grouped matrix-valued time series, separating global and local main effects from bilinear common components.
