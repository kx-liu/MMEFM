.standardize_detection_residuals <- function(check_Y) {
  scale <- numeric(length(check_Y))
  names(scale) <- names(check_Y)
  check_Z <- check_Y
  for (m in seq_along(check_Y)) {
    scale[m] <- sqrt(mean(check_Y[[m]]^2))
    if (!is.finite(scale[m]) || scale[m] <= 0) {
      group <- if (is.null(names(check_Y))) m else names(check_Y)[m]
      stop("Residual RMS scale must be finite and strictly positive in group ", group, ".")
    }
    check_Z[[m]] <- check_Y[[m]] / scale[m]
  }
  list(check_Z = check_Z, scale = scale)
}

.global_detection_statistic <- function(check_Z, rank, Q_hat, J_hat) {
  M <- length(check_Z)
  TT <- dim(check_Z[[1L]])[1L]
  Z <- vector("list", M)
  denominator <- numeric(M)
  for (m in seq_len(M)) {
    p <- dim(check_Z[[m]])[2L]
    q <- dim(check_Z[[m]])[3L]
    row_cov <- matrix(0, p, p)
    col_cov <- matrix(0, q, q)
    Z[[m]] <- matrix(0, TT, as.numeric(rank$kr) * rank$kc)
    for (t in seq_len(TT)) {
      Yt <- matrix(check_Z[[m]][t, , , drop = FALSE], p, q)
      Z[[m]][t, ] <- as.vector(crossprod(Q_hat[[m]], Yt) %*% J_hat[[m]])
      row_cov <- row_cov + tcrossprod(Yt) / TT
      col_cov <- col_cov + crossprod(Yt) / TT
    }
    U_hat <- .leading_eigenvectors(row_cov, rank$kr + rank$kr_m[m])
    V_hat <- .leading_eigenvectors(col_cov, rank$kc + rank$kc_m[m])
    energy <- matrix(0, ncol(U_hat), ncol(V_hat))
    for (t in seq_len(TT)) {
      Yt <- matrix(check_Z[[m]][t, , , drop = FALSE], p, q)
      energy <- energy + (crossprod(U_hat, Yt) %*% V_hat)^2 / TT
    }
    denominator[m] <- max(energy)
    if (!is.finite(denominator[m]) || denominator[m] <= 0) {
      group <- if (is.null(names(check_Z))) m else names(check_Z)[m]
      stop("Detection denominator must be finite and strictly positive in group ", group, ".")
    }
  }
  SG <- numeric(M)
  names(SG) <- names(check_Z)
  for (m in seq_len(M)) {
    numerator <- 0
    for (n in seq_len(M)) {
      if (n != m) numerator <- numerator + max((crossprod(Z[[m]], Z[[n]]) / TT)^2)
    }
    SG[m] <- numerator / denominator[m]^2
    if (!is.finite(SG[m])) stop("Nonfinite detection statistic in group ", if (is.null(names(check_Z))) m else names(check_Z)[m], ".")
  }
  list(statistic = SG, second_largest = unname(sort(SG, decreasing = TRUE)[2L]))
}

.circular_shift_detection_groups <- function(check_Z, reference_group) {
  TT <- dim(check_Z[[1L]])[1L]
  shifted <- check_Z
  for (m in seq_along(check_Z)) {
    if (m != reference_group) {
      shift <- sample.int(TT, 1L) - 1L
      index <- 1L + ((seq_len(TT) - 1L + as.numeric(shift)) %% TT)
      shifted[[m]] <- check_Z[[m]][index, , , drop = FALSE]
    }
  }
  shifted
}

.bootstrap_global_detection <- function(check_Z, rank, B, reference_group, K0, max_iter, tol, verbose) {
  second_largest <- numeric(B)
  nonconverged <- 0L
  for (b in seq_len(B)) {
    draw <- tryCatch({
      shifted <- .circular_shift_detection_groups(check_Z, reference_group)
      loadings <- .initial_common_loadings(shifted, rank$kr, rank$kc, K0, max_iter, tol, verbose)
      statistic <- .global_detection_statistic(shifted, rank, loadings$Q_hat, loadings$J_hat)
      list(second_largest = statistic$second_largest, converged = loadings$converged)
    }, error = function(e) stop("Bootstrap replication ", b, " failed: ", conditionMessage(e), call. = FALSE))
    second_largest[b] <- draw$second_largest
    if (!draw$converged) nonconverged <- nonconverged + 1L
  }
  list(second_largest = second_largest, nonconverged = nonconverged)
}

#' Detect and screen global common factors in MMEFM
#'
#' Practical manuscript detection based on group statistics S_G,m and circular-shift bootstrap calibration.
#' @param Xt A canonical list of finite T x p_m x q_m arrays (see [est_MMEFM()]), or an `mmefm_fit`. Raw arrays are reduced to the direct main-effect residual `check_Y`; a fit supplies its stored `check_Y`, whose dimensions are validated anew.
#' @param rank NULL to select ranks for raw data or use a supplied fit's ranks; otherwise a complete canonical eight-field rank list or `mmefm_rank`. Supplied ranks override fit ranks. Statistical values are validated again; rank-object diagnostics are preserved.
#' @param B Positive whole-number integer-range bootstrap size, default 199. Re-estimation is expensive; use small values such as 3 for examples and smoke checks.
#' @param alpha Finite calibration level strictly between zero and one.
#' @param reference_group Integer group index between 1 and M, held fixed during bootstrap.
#' @inheritParams select_MMEFM_rank
#' @details Each residual group is divided by its RMS scale, `sqrt(mean(check_Y^2))`, without further centering. Algorithm 1 re-estimates initial Q/J spaces on these standardized arrays, even for fit input; stored fitted loadings are never reused.
#'
#' S_G,m sums, over other groups, the largest squared cross covariance between projected global-core coordinates, divided by the square of the strongest within-group projected-coordinate mean square. Denominator row/column spaces use total global plus local common ranks, ignoring their split.
#'
#' Algorithm 2 holds the reference group fixed and independently circularly shifts every other group's whole time series uniformly by 0 through T-1. Every replication re-estimates Q/J and denominator loading spaces at fixed ranks. One Type-8 empirical (1-alpha) quantile of bootstrap second-largest S_G values calibrates both existence and screening. Existence requires the observed second-largest statistic to strictly exceed the cutoff. Final screening is empty unless existence rejects and at least two groups strictly exceed that same cutoff.
#'
#' All detector randomness uses one scoped seed after any separately scoped rank selection; caller RNG state is restored on success and error. Observed projection non-convergence warns once; bootstrap non-convergence warns once in aggregate. Final iterates remain usable. Invalid scales, denominators, and numerical failures error without floors, repair, or discarded draws; bootstrap errors identify the replication. Whole-sample circular shifts assume a single stationary segment.
#' @return An `mmefm_global_detection` list containing exactly `call`, group `statistic`, observed `second_largest`, `cutoff`, logical `exists`, group logical vectors `initial_selected` and `selected`, group RMS `scale`, canonical `rank`, and `bootstrap`. The bootstrap list contains exactly `B`, `alpha`, `reference_group`, and the length-B numeric calibration sample `second_largest`. Group-indexed vectors preserve optional names. Loading matrices and residual workspaces are not returned.
#' @export
detect_MMEFM_global <- function(Xt, rank = NULL, B = 199L, alpha = 0.05, reference_group = 1L,
                               K0 = 20L, max_iter = 20L, tol = 1e-4, seed = 2026L, verbose = FALSE) {
  call <- match.call()
  fit_input <- inherits(Xt, "mmefm_fit")
  data_info <- .validate_Xt(if (fit_input) Xt$check_Y else Xt)
  controls <- list(B = B, reference_group = reference_group, K0 = K0, max_iter = max_iter, seed = seed)
  for (name in names(controls)) {
    value <- controls[[name]]
    minimum <- if (name == "seed") 0L else 1L
    maximum <- if (name == "reference_group") min(data_info$M, .Machine$integer.max) else .Machine$integer.max
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != floor(value) || value < minimum || value > maximum) {
      stop(name, " must be one finite whole number between ", minimum, " and ", maximum, ".")
    }
  }
  B <- as.integer(B)
  reference_group <- as.integer(reference_group)
  K0 <- as.integer(K0)
  max_iter <- as.integer(max_iter)
  seed <- as.integer(seed)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be one finite number strictly between zero and one.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be one finite positive number.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be one nonmissing logical.")
  if (fit_input) {
    check_Y <- Xt$check_Y
    if (is.null(rank)) rank <- Xt$rank
  } else {
    check_Y <- .main_effect_moments(Xt)$check_Y
    if (is.null(rank)) rank <- select_MMEFM_rank(Xt, K0, max_iter, tol, seed, verbose)
  }
  ranks <- .validate_rank(rank, data_info)
  rank <- .new_mmefm_rank(ranks, if (inherits(rank, "mmefm_rank")) rank$diagnostics else NULL)
  standardized <- .standardize_detection_residuals(check_Y)
  computed <- .with_preserved_seed(seed, {
    loadings <- .initial_common_loadings(standardized$check_Z, ranks$kr, ranks$kc, K0, max_iter, tol, verbose)
    if (!loadings$converged) warning("Observed detection projection reached max_iter without convergence.", call. = FALSE)
    observed <- .global_detection_statistic(standardized$check_Z, ranks, loadings$Q_hat, loadings$J_hat)
    bootstrap <- .bootstrap_global_detection(standardized$check_Z, ranks, B, reference_group, K0, max_iter, tol, verbose)
    if (bootstrap$nonconverged > 0L) warning(bootstrap$nonconverged, " of ", B, " bootstrap projections reached max_iter without convergence.", call. = FALSE)
    list(observed = observed, bootstrap = bootstrap)
  })
  cutoff <- stats::quantile(computed$bootstrap$second_largest, 1 - alpha, names = FALSE, type = 8L)
  observed <- computed$observed
  exists <- observed$second_largest > cutoff
  initial_selected <- observed$statistic > cutoff
  selected <- initial_selected
  if (!exists || sum(initial_selected) < 2L) selected[] <- FALSE
  structure(list(call = call, statistic = observed$statistic, second_largest = observed$second_largest,
                 cutoff = cutoff, exists = exists, initial_selected = initial_selected, selected = selected,
                 scale = standardized$scale, rank = rank,
                 bootstrap = list(B = B, alpha = alpha, reference_group = reference_group, second_largest = computed$bootstrap$second_largest)),
            class = "mmefm_global_detection")
}

#' Print global common-factor detection
#' @param x An `mmefm_global_detection` object.
#' @param ... Further arguments, currently unused.
#' @return The input invisibly.
#' @export
print.mmefm_global_detection <- function(x, ...) {
  cat(sprintf("MMEFM global detection: %d groups\n", length(x$statistic)))
  cat(sprintf("Second-largest S_G: %.6g; cutoff: %.6g; alpha: %.3g\n", x$second_largest, x$cutoff, x$bootstrap$alpha))
  cat("Global common factors detected:", if (x$exists) "yes" else "no", "\n")
  groups <- which(x$selected)
  labels <- if (is.null(names(x$selected))) as.character(groups) else names(x$selected)[groups]
  cat("Selected groups:", if (length(labels)) paste(labels, collapse = ", ") else "none", "\n")
  invisible(x)
}
