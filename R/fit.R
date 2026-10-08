#' Estimate a Multilevel Main Effects Matrix Factor Model
#'
#' Fits global and local main effects and common matrix-factor components at complete positive ranks.
#' @inheritParams select_MMEFM_rank
#' @param rank NULL for automatic rank selection, a complete named list with `r1`, `l1`, `r2`, `l2`, `kr`, `kc`, `kr_m`, and `kc_m`, or an `mmefm_rank` object. Global ranks are positive scalar whole numbers; local ranks are positive whole-number vectors in group order. Named local ranks are reordered to named `Xt` groups. Each rank is smaller than its corresponding spatial dimension, and under IC1, `r1 + r2[m]` and `kr + kr_m[m]` must be at most `p_m - 1`, while `l1 + l2[m]` and `kc + kc_m[m]` must be at most `q_m - 1`. Equality at dimension minus one is allowed.
#' @param alignment_method Global common-coordinate alignment: `"VanLoan"` (default) or `"Procrustes"`.
#' @param max_iter_procrustes Positive whole-number integer-range limit for Procrustes alignment iterations.
#' @param tol_procrustes Finite positive relative-objective tolerance for Procrustes alignment.
#' @param refit_method Global common-loading refit: `"VanLoan"` (default) or `"ALS"`. Alignment and refit choices are independent.
#' @param max_iter_als Positive whole-number integer-range limit for ALS refit iterations.
#' @param tol_als Finite positive relative-residual-objective tolerance for ALS refitting.
#' @param lambda Finite nonnegative ridge parameter for main-effect Gram systems only. The default zero implements the manuscript estimator without regularization.
#' @details Automatic selection is followed by a fresh known-rank fit with the same seed, so supplying the selected rank explicitly gives the same statistical estimates. Both public stochastic stages preserve the caller's RNG state, including on error.
#'
#' Global main-effect ranks must not exceed T. If either alignment or refitting uses Van Loan, `kr * kc <= T` is also necessary. These dimension checks do not guarantee nonsingular Gram systems: singular unregularized systems remain visible errors, with no numerical fallback or silent rank repair.
#'
#' Iteration-limit non-convergence produces a warning for each affected active stage: rank-selection projection, final projection, Procrustes alignment, or ALS refitting. Final iterates are retained, with small operational convergence metadata available programmatically.
#' @return An `mmefm_fit` list with exactly `call`, `dimensions`, `rank` (an `mmefm_rank`), `main_effect`, `common_component`, `check_Y`, and `convergence`. Observed arrays are not stored separately.
#'
#' Main effects retain `direct`, `global`, and `local` estimates. Common components retain `global`, `local`, and method-specific `alignment` estimates. Global `Q_hat`/`J_hat` are the original Algorithm 1 loadings, while `Q_tilde`/`J_tilde` are final refitted loadings. `G_hat` is estimated after alignment; `G_tilde` equals `G_hat` for Van Loan refitting and is updated after every Q/J sweep for ALS. Local common estimates are `R_hat`, `C_hat`, and `F_hat`.
#'
#' `convergence` contains final-estimation `projection`, `procrustes`, and `als` iteration counts and logical indicators; inactive iterative methods use NULL. See [fitted.mmefm_fit()] for component reconstruction and [residuals.mmefm_fit()] for residuals.
#' @export
est_MMEFM <- function(
    Xt, rank = NULL, K0 = 20L, max_iter = 20L, tol = 1e-4,
    alignment_method = c("VanLoan", "Procrustes"), max_iter_procrustes = 20L, tol_procrustes = 1e-4,
    refit_method = c("VanLoan", "ALS"), max_iter_als = 20L, tol_als = 1e-4,
    lambda = 0, seed = 2026L, verbose = FALSE) {
  call <- match.call()
  data_info <- .validate_Xt(Xt)
  alignment_method <- match.arg(alignment_method)
  refit_method <- match.arg(refit_method)
  integer_controls <- list(K0 = K0, max_iter = max_iter, max_iter_procrustes = max_iter_procrustes,
                           max_iter_als = max_iter_als, seed = seed)
  for (name in names(integer_controls)) {
    value <- integer_controls[[name]]
    minimum <- if (name == "seed") 0L else 1L
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum || value > .Machine$integer.max) {
      stop(sprintf("%s must be one finite whole number >= %d within R's integer range.", name, minimum))
    }
  }
  K0 <- as.integer(K0)
  max_iter <- as.integer(max_iter)
  max_iter_procrustes <- as.integer(max_iter_procrustes)
  max_iter_als <- as.integer(max_iter_als)
  seed <- as.integer(seed)
  for (name in c("tol", "tol_procrustes", "tol_als")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(sprintf("%s must be one finite positive number.", name))
    }
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) stop("lambda must be one finite nonnegative number.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be one nonmissing logical.")
  if (is.null(rank)) rank <- select_MMEFM_rank(Xt, K0, max_iter, tol, seed, verbose)
  ranks <- .validate_rank(rank, data_info)
  diagnostics <- if (inherits(rank, "mmefm_rank")) rank$diagnostics else NULL
  rank <- .new_mmefm_rank(ranks, diagnostics)
  if (ranks$r1 > data_info$T || ranks$l1 > data_info$T) stop("Global main-effect ranks r1 and l1 must not exceed T.")
  if ((alignment_method == "VanLoan" || refit_method == "VanLoan") &&
      as.numeric(ranks$kr) * ranks$kc > data_info$T) {
    stop("VanLoan alignment or refitting requires kr * kc <= T.")
  }
  main <- .estimate_main_effect(Xt, ranks, lambda)
  common <- .with_preserved_seed(seed,
    .estimate_common_component(main$check_Y, ranks, K0, max_iter, tol,
                               alignment_method, max_iter_procrustes, tol_procrustes,
                               refit_method, max_iter_als, tol_als, verbose))
  convergence <- common$convergence
  if (!convergence$projection$converged) warning("Common-loading projection reached max_iter without convergence.", call. = FALSE)
  if (alignment_method == "Procrustes") {
    failed <- which(!convergence$procrustes$converged)
    failed <- setdiff(failed, common$common_component$alignment$reference_group)
    if (length(failed)) {
      groups <- if (is.null(data_info$group_names)) failed else data_info$group_names[failed]
      warning("Procrustes alignment reached max_iter_procrustes without convergence for group(s): ",
              paste(groups, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (refit_method == "ALS" && !convergence$als$converged) warning("ALS refitting reached max_iter_als without convergence.", call. = FALSE)
  structure(list(call = call, dimensions = data_info, rank = rank, main_effect = main$main_effect,
                 common_component = common$common_component, check_Y = main$check_Y, convergence = convergence),
            class = "mmefm_fit")
}

#' Reconstruct fitted MMEFM components
#'
#' @param object An `mmefm_fit` object.
#' @param component Component to reconstruct: `"global_main"` uses the global mean and refitted main loadings; `"local_main"` uses local means, loadings, and scores; `"main"` sums these structured main effects; `"global_common"` uses `Q_tilde`, `G_tilde`, and `J_tilde`; `"local_common"` uses `R_hat`, `F_hat`, and `C_hat`; `"all"` sums main and both common components.
#' @param ... Further arguments, currently unused.
#' @return A list of T x p_m x q_m arrays preserving group names. Main components use the final structured estimator, not the raw direct moments. `check_Y` is not an additional fitted component.
#' @export
fitted.mmefm_fit <- function(object, component = c("all", "main", "global_main", "local_main", "global_common", "local_common"), ...) {
  component <- match.arg(component)
  dimensions <- object$dimensions
  result <- vector("list", dimensions$M)
  names(result) <- dimensions$group_names
  global_main <- object$main_effect$global
  local_main <- object$main_effect$local
  global_common <- object$common_component$global
  local_common <- object$common_component$local
  for (m in seq_len(dimensions$M)) {
    p <- unname(dimensions$p[m])
    q <- unname(dimensions$q[m])
    one_p <- matrix(1, p, 1L)
    one_q <- matrix(1, q, 1L)
    result[[m]] <- array(0, c(dimensions$T, p, q))
    for (t in seq_len(dimensions$T)) {
      value <- matrix(0, p, q)
      if (component %in% c("all", "main", "global_main")) {
        row_effect <- global_main$A1_tilde[[m]] %*% t(global_main$alpha_hat[t, , drop = FALSE])
        column_effect <- global_main$B1_tilde[[m]] %*% t(global_main$beta_hat[t, , drop = FALSE])
        value <- value + global_main$mu_hat[t, 1L] * tcrossprod(one_p, one_q) +
          tcrossprod(row_effect, one_q) + tcrossprod(one_p, column_effect)
      }
      if (component %in% c("all", "main", "local_main")) {
        row_effect <- local_main$A2_hat[[m]] %*% t(local_main$alpha_hat[[m]][t, , drop = FALSE])
        column_effect <- local_main$B2_hat[[m]] %*% t(local_main$beta_hat[[m]][t, , drop = FALSE])
        value <- value + local_main$mu_hat[[m]][t, 1L] * tcrossprod(one_p, one_q) +
          tcrossprod(row_effect, one_q) + tcrossprod(one_p, column_effect)
      }
      if (component %in% c("all", "global_common")) {
        Gt <- matrix(global_common$G_tilde[t, , , drop = FALSE], object$rank$kr, object$rank$kc)
        value <- value + global_common$Q_tilde[[m]] %*% Gt %*% t(global_common$J_tilde[[m]])
      }
      if (component %in% c("all", "local_common")) {
        Ft <- matrix(local_common$F_hat[[m]][t, , , drop = FALSE], object$rank$kr_m[m], object$rank$kc_m[m])
        value <- value + local_common$R_hat[[m]] %*% Ft %*% t(local_common$C_hat[[m]])
      }
      result[[m]][t, , ] <- value
    }
  }
  result
}

#' Residuals of an MMEFM fit
#'
#' @inheritParams fitted.mmefm_fit
#' @return A named list of T x p_m x q_m residual arrays. Observations are reconstructed from direct `c_hat`, `a_hat`, `b_hat`, and `check_Y`, then the complete final structured fit is subtracted. No separate observed-data copy is stored.
#' @export
residuals.mmefm_fit <- function(object, ...) {
  result <- object$check_Y
  direct <- object$main_effect$direct
  fitted_values <- fitted.mmefm_fit(object, "all")
  for (m in seq_len(object$dimensions$M)) {
    p <- object$dimensions$p[m]
    q <- object$dimensions$q[m]
    for (t in seq_len(object$dimensions$T)) {
      observed <- matrix(result[[m]][t, , , drop = FALSE], p, q) + matrix(direct$c_hat[[m]][t, 1L], p, q)
      observed <- sweep(observed, 1L, direct$a_hat[[m]][t, ], "+")
      observed <- sweep(observed, 2L, direct$b_hat[[m]][t, ], "+")
      result[[m]][t, , ] <- observed - matrix(fitted_values[[m]][t, , , drop = FALSE], p, q)
    }
  }
  result
}

#' Summarize an MMEFM fit
#'
#' @inheritParams fitted.mmefm_fit
#' @return A compact object of class `summary.mmefm_fit` containing the call, dimensions, eight ranks without diagnostics, alignment method, and refit method. Operational convergence metadata is omitted.
#' @export
summary.mmefm_fit <- function(object, ...) {
  structure(list(call = object$call, dimensions = object$dimensions,
                 rank = unclass(object$rank)[c("r1", "l1", "r2", "l2", "kr", "kc", "kr_m", "kc_m")],
                 alignment_method = object$common_component$alignment$method,
                 refit_method = if (is.null(object$convergence$als)) "VanLoan" else "ALS"),
            class = "summary.mmefm_fit")
}

#' @rdname summary.mmefm_fit
#' @param x An `mmefm_fit` or `summary.mmefm_fit` object, as appropriate.
#' @return Print methods return their input invisibly.
#' @export
print.summary.mmefm_fit <- function(x, ...) {
  cat(sprintf("MMEFM fit: %d groups, T = %d\n", x$dimensions$M, x$dimensions$T))
  global <- x$rank[c("r1", "l1", "kr", "kc")]
  cat("Global ranks:", paste(paste(names(global), unlist(global), sep = "="), collapse = ", "), "\n")
  cat(sprintf("Alignment: %s; refit: %s\n", x$alignment_method, x$refit_method))
  local <- data.frame(p = x$dimensions$p, q = x$dimensions$q,
                      r2 = x$rank$r2, l2 = x$rank$l2, kr_m = x$rank$kr_m, kc_m = x$rank$kc_m)
  rownames(local) <- if (is.null(x$dimensions$group_names)) as.character(seq_len(x$dimensions$M)) else x$dimensions$group_names
  print(local)
  invisible(x)
}

#' @rdname summary.mmefm_fit
#' @export
print.mmefm_fit <- function(x, ...) {
  print.summary.mmefm_fit(summary.mmefm_fit(x))
  invisible(x)
}
