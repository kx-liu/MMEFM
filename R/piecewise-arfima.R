#' Generate piecewise ARFIMA time series
#'
#' Simulate independent ARFIMA(p, d, 0) series with a stable AR short-memory component and segment-specific fractional memory. The memory specification may be shared across series or differ by series.
#' @param T Positive whole-number sample length within R's integer range.
#' @param n Positive whole-number number of series within R's integer range. The series dimension is retained when n = 1.
#' @param d_vec A shared numeric vector of finite memory parameters below 0.5, a shared `list(d = ..., seg_lens = ...)`, or a list of length n specifying each series. Numeric vectors split T approximately equally, with the remainder in the final segment. Explicit segment lengths are positive whole numbers summing to T. No lower bound on d is imposed.
#' @param phi Nonempty finite AR coefficient vector in lag order: `X[t] = phi[1]*X[t-1] + ... + phi[p]*X[t-p] + innovation[t]`. All roots of `1 - phi[1]*z - ... - phi[p]*z^p` must lie strictly outside the unit circle. A scalar specifies AR(1) with absolute value below one; stable higher-order coefficients may individually exceed one in magnitude.
#' @param innovation Gaussian (default) or raw, unstandardized Student-t innovations.
#' @param df Finite positive Student-t degrees of freedom, used for Student-t draws only. Small df may give infinite variance or violate the manuscript's moment assumptions.
#' @param burn Nonnegative whole-number burn-in within R's integer range.
#' @param filter_length NULL for `min(2000, max(200, 5*T))`, or a nonnegative whole-number filter truncation length within R's integer range.
#' @param seed Nonnegative whole-number seed within R's integer range. The caller's RNG state and kind are preserved on success and error; repeated calls agree under the same seed and RNG kind.
#' @param return_filters One nonmissing logical. TRUE returns deterministic filter information alongside the series.
#' @details
#' Within a segment of memory d, the fractional filter represents \eqn{(1-B)^{-d}}, where B is the lag operator. Zero memory gives an AR process; positive memory below 0.5 introduces long-memory dependence. Changing d between segments need not give a stationary series over the whole sample.
#'
#' Nonzero memory uses a truncated AR impulse response followed by segment-specific fractional filters. The output segment determines its fractional parameter; the innovation stream is not restarted at a boundary. Entirely zero-memory series use direct AR recursion from zero pre-sample values. Scalar AR(1) uses the corresponding geometric filter or recursion. The retained sample follows `burn + 2*filter_length` warmup observations, and innovations are generated in series-column order.
#'
#' AR stability concerns only the short-memory component. The bounds on d and df do not by themselves establish every theoretical assumption. Nonfinite generated series error without numerical repair.
#' @return With `return_filters = FALSE`, a finite numeric `T x n` matrix. With TRUE, a list containing `series` (the same matrix), supplied `phi`, `filter_length`, `short_filter`, and `series_specs`.
#'
#' `short_filter` contains AR impulse-response coefficients at lags 0 through `filter_length`. Each of the n `series_specs` contains `d`, `seg_lens`, inclusive `seg_start`/`seg_end` indices, and one `fractional_filters` coefficient vector per segment at those same lags. For zero memory, `short_filter` is an inspection representation; simulation uses recursion. Returning metadata consumes no additional random draws.
#' @seealso [gen_MMEFM()] for grouped matrix-valued simulation
#' @examples
#' series <- gen_piecewise_arfima(
#'   T = 12L, n = 2L, d_vec = c(0, 0.2), phi = c(0.4, -0.1),
#'   burn = 4L, filter_length = 8L, seed = 2026L
#' )
#' dim(series)
#' @export
gen_piecewise_arfima <- function(
    T, n = 1L, d_vec = 0, phi = 0.5, innovation = c("Gaussian", "Student-t"),
    df = 6, burn = 200L, filter_length = NULL, seed = 2026L, return_filters = FALSE) {
  innovation <- match.arg(innovation)
  controls <- list(T = T, n = n, burn = burn, seed = seed)
  if (!is.null(filter_length)) {
    controls$filter_length <- filter_length
  }
  for (name in names(controls)) {
    value <- controls[[name]]
    minimum <- if (name %in% c("T", "n")) 1L else 0L
    if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum || value > .Machine$integer.max) {
      stop(name, " must be one finite whole number within integer range, >= ", minimum, ".")
    }
  }
  T <- as.integer(T)
  n <- as.integer(n)
  burn <- as.integer(burn)
  seed <- as.integer(seed)
  if (is.null(filter_length)) {
    filter_length <- min(2000, max(200, 5 * as.numeric(T)))
  }
  filter_length <- as.integer(filter_length)
  .validate_simulation_ar(phi, "phi")
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("df must be one finite positive number.")
  }
  if (!is.logical(return_filters) || length(return_filters) != 1L || is.na(return_filters)) {
    stop("return_filters must be one nonmissing logical.")
  }

  shared <- is.numeric(d_vec) ||
    (is.list(d_vec) && !is.null(names(d_vec)) && any(names(d_vec) %in% c("d", "seg_lens")))
  if (shared) {
    spec <- .simulation_memory_spec(d_vec, T, "d_vec")
    specs <- rep(list(spec), n)
  } else {
    if (!is.list(d_vec) || is.data.frame(d_vec) || length(d_vec) != n) {
      stop("d_vec must be a shared memory specification or a list of n series specifications.")
    }
    specs <- vector("list", n)
    for (j in seq_len(n)) {
      specs[[j]] <- .simulation_memory_spec(d_vec[[j]], T, paste0("d_vec[[", j, "]]"))
    }
  }

  series <- .with_preserved_seed(seed, {
    if (shared) {
      .simulation_series(n, T, spec, phi, innovation, df, burn, filter_length)
    } else {
      result <- matrix(0, T, n)
      for (j in seq_len(n)) {
        result[, j] <- .simulation_series(1L, T, specs[[j]], phi, innovation, df, burn, filter_length)
      }
      result
    }
  })
  if (any(!is.finite(series))) {
    stop("Nonfinite generated piecewise ARFIMA series.")
  }
  if (!return_filters) {
    return(series)
  }

  # The same filter recurrences serve simulation and deterministic inspection.
  for (j in seq_len(n)) {
    end <- cumsum(specs[[j]]$seg_lens)
    specs[[j]]$seg_start <- c(1L, end[-length(end)] + 1L)
    specs[[j]]$seg_end <- end
    filters <- vector("list", length(specs[[j]]$d))
    for (segment in seq_along(filters)) {
      filters[[segment]] <- .simulation_fractional_filter(specs[[j]]$d[segment], filter_length)
    }
    specs[[j]]$fractional_filters <- filters
  }
  list(
    series = series, phi = phi, filter_length = filter_length,
    short_filter = .simulation_ar_filter(phi, filter_length), series_specs = specs
  )
}
