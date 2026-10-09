## Release status

First submission of MMEFM, version 0.1.0. This is a release candidate; the package has not been submitted to CRAN.

## Test environments

* Local: macOS 27.0.1, R 4.6.1 (2026-06-24), aarch64-apple-darwin23; full `R CMD check --as-cran`.
* GitHub Actions: macOS 26.6.2 (aarch64), Windows x64 (build 26100), and Ubuntu 24.04.5 LTS (x86_64), each with R release 4.6.1; Ubuntu 24.04.5 LTS with R-devel (2026-10-08 r90650). These checks use `R CMD check --no-manual`, not `--as-cran`.

The four-platform results were verified from [workflow run 37982998797](https://github.com/kx-liu/MMEFM/actions/runs/37982998797), targeting code commit `45f0c8ac545b62de64fff7440b14e9a9cf3d10e2`.

## R CMD check results

`R CMD check --as-cran` on `MMEFM_0.1.0.tar.gz`: 0 ERRORs, 0 WARNINGs, 1 NOTE.

* CRAN incoming feasibility: "New submission". This is expected for the first submission; no development-version NOTE remains.

HTML validation, including mathematical rendering, completed successfully with V8 8.2.0 and HTML Tidy 5.8.0. Evaluated examples, vignette checks and rebuilding, and the PDF manual passed. All 2,316 regression-test assertions passed, with 0 failures, warnings, or skips.

All four GitHub Actions checks reported `Status: OK`, with 0 ERRORs, 0 WARNINGs, and 0 NOTEs. Each passed the examples, vignette checks, and all 2,316 regression-test assertions, with 0 failures, warnings, or skips.

## Additional remarks

MMEFM implements statistical methodology for grouped matrix-valued time series, separating global and local main effects from bilinear common components.
