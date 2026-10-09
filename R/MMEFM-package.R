#' Multilevel Main Effects Matrix Factor Model Tools
#'
#' Models grouped matrix-valued time series by distinguishing global and local dynamics across groups and separating grand means and row and column main effects from bilinear interaction components.
#' @details [gen_MMEFM()] simulates grouped data and generating quantities; [gen_piecewise_arfima()] simulates scalar series with segment-specific memory. [select_MMEFM_rank()] estimates loading ranks, and [est_MMEFM()] fits the model at supplied or automatically selected ranks. Use `summary()`, `fitted()`, and `residuals()` to inspect a fit and reconstruct its components. [detect_MMEFM_global()] tests for a global common component and screens the groups sharing it.
#'
#' The accompanying methodology is developed in *Multilevel Main Effects Matrix Factor Model* by Kaixin Liu, Clifford Lam, and Zetai Cen. Main effects may be nonstationary under the manuscript's conditions. Whole-sample circular-shift detection, however, requires a single stationary segment. See `vignette("getting-started", package = "MMEFM")` for a worked example.
#' @keywords internal
"_PACKAGE"
