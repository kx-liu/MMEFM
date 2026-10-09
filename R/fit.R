#' Estimate a Multilevel Main Effects Matrix Factor Model
#'
#' Estimate global and local grand means, row and column main effects, and common matrix-factor components in grouped matrix-valued time series.
#' @inheritParams select_MMEFM_rank
#' @param rank NULL for automatic selection, a complete named list containing exactly `r1`, `l1`, `r2`, `l2`, `kr`, `kc`, `kr_m`, and `kc_m`, or an `mmefm_rank` object. Global ranks are positive scalar whole numbers; local ranks are positive whole-number vectors in group order. Named local vectors are matched to named `Xt` groups. Each rank is smaller than its spatial dimension. Under IC1, row sums `r1 + r2[m]` and `kr + kr_m[m]` are at most `p_m - 1`; column sums `l1 + l2[m]` and `kc + kc_m[m]` are at most `q_m - 1`. Equality is allowed.
#' @param alignment_method Method for aligning group-specific global common-factor coordinates: `"VanLoan"` (default) or `"Procrustes"`.
#' @param max_iter_procrustes Positive whole-number limit on Procrustes alignment iterations.
#' @param tol_procrustes Finite positive tolerance for relative change in the Procrustes alignment score.
#' @param refit_method Method for refitting global common loadings: `"VanLoan"` (default) or `"ALS"`. This choice is independent of `alignment_method`.
#' @param max_iter_als Positive whole-number limit on alternating least-squares (ALS) refitting iterations.
#' @param tol_als Finite positive tolerance for relative change in the squared-residual objective during ALS refitting.
#' @param lambda Finite nonnegative ridge parameter for main-effect Gram systems only. Zero, the default, gives the unregularized manuscript estimator.
#' @details
#' Direct moment estimators first recover the combined grand mean and centered row and column main effects in each group. Cross-group moments and projections separate global and local main-effect structures. Shared main-effect scores are estimated by pooling temporal information, and their global loadings are refitted.
#'
#' The direct additive estimates are subtracted to form `check_Y`. Iterative projection estimates initial global common loading spaces from these residuals. Projections then recover local common loadings and factors. After removing local common components, group-specific global common-factor coordinates are aligned and pooled into one shared factor series. Global common loadings are refitted to that series.
#'
#' Van Loan alignment approximates the between-group coordinate mapping by a Kronecker product. Van Loan refitting similarly approximates the unconstrained fitted loading product, with the aligned global core held fixed. Optional Procrustes alignment alternates orthogonal row and column transformations. Optional ALS refitting alternates row and column loading updates and updates the shared global core after each sweep.
#'
#' The direct main-effect estimates are not the final low-rank main-effect reconstruction. [fitted.mmefm_fit()] uses global and local factor scores and loadings, including refitted global main loadings. `check_Y` is an intermediate residual, not the final model residual; use [residuals.mmefm_fit()] for the latter.
#'
#' With `rank = NULL`, [select_MMEFM_rank()] is followed by a fresh known-rank fit using the same seed. Supplying that selected rank explicitly gives the same statistical estimates. Both stochastic stages preserve the caller's RNG state, including on error.
#'
#' Global main-effect ranks must not exceed T. Van Loan alignment or refitting additionally requires `kr * kc <= T`. These bounds do not guarantee nonsingular Gram systems. Singular unregularized systems error without numerical fallback or rank repair. Iteration-limit warnings retain the final iterates; inspect `fit$convergence`.
#' @return An object of class `mmefm_fit` with `call`, `dimensions`, `rank`, `main_effect`, `common_component`, `check_Y`, and `convergence`. `rank` is an `mmefm_rank`; observed arrays are not stored separately.
#'
#' `main_effect` contains `direct` estimates (`c_hat`, `a_hat`, `b_hat`), plus `global` and `local` grand means, loadings, and factor scores. Initial global loadings `A1_hat` and `B1_hat` are retained alongside refitted `A1_tilde` and `B1_tilde`. Local loadings are `A2_hat` and `B2_hat`; both levels store `alpha_hat` and `beta_hat`.
#'
#' `common_component` contains `global`, `local`, and method-specific `alignment` estimates. Global `Q_hat` and `J_hat` are initial loading estimates; `Q_tilde` and `J_tilde` are refitted loadings. `G_hat` is the aligned, pooled core. `G_tilde` equals `G_hat` for Van Loan refitting and is updated during ALS. Local estimates are `R_hat`, `C_hat`, and `F_hat`.
#'
#' `check_Y` contains the direct additive residual arrays. `convergence` records iteration counts and convergence indicators for final `projection`, `procrustes`, and `als` stages; inactive methods are NULL. [summary.mmefm_fit()] reports dimensions, ranks, and method choices, not convergence diagnostics.
#' @seealso [gen_MMEFM()], [select_MMEFM_rank()], [detect_MMEFM_global()], [fitted.mmefm_fit()], [residuals.mmefm_fit()]
#' @examples
#' rank <- list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3), l2 = rep(1L, 3),
#'              kr = 1L, kc = 1L, kr_m = rep(1L, 3), kc_m = rep(1L, 3))
#' simulation <- gen_MMEFM(T = 40L, p = c(a = 6L, b = 7L, c = 8L), q = c(a = 7L, b = 8L, c = 6L),
#'                         rank = rank, burn = 20L, filter_length = 40L, seed = 2026L)
#' fit <- est_MMEFM(simulation$Xt, rank = simulation$rank, K0 = 8L, max_iter = 50L, seed = 2026L)
#' summary(fit)
#' fit$convergence
#' lapply(fitted(fit), dim)
#' lapply(fitted(fit, component = "global_common"), dim)
#' vapply(residuals(fit), function(x) sqrt(mean(x^2)), numeric(1L))
#' @export
est_MMEFM <- function(
    Xt, rank = NULL, alignment_method = c("VanLoan", "Procrustes"), refit_method = c("VanLoan", "ALS"),
    K0 = 20L, max_iter = 20L, tol = 1e-4, max_iter_procrustes = 20L, tol_procrustes = 1e-4,
    max_iter_als = 20L, tol_als = 1e-4,
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
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite nonnegative number.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one nonmissing logical.")
  }
  if (is.null(rank)) {
    rank <- select_MMEFM_rank(Xt, K0 = K0, max_iter = max_iter, tol = tol, seed = seed, verbose = verbose)
  }
  ranks <- .validate_rank(rank, data_info)
  diagnostics <- if (inherits(rank, "mmefm_rank")) rank$diagnostics else NULL
  rank <- .new_mmefm_rank(ranks, diagnostics)
  if (ranks$r1 > data_info$T || ranks$l1 > data_info$T) {
    stop("Global main-effect ranks r1 and l1 must not exceed T.")
  }
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
  if (!convergence$projection$converged) {
    warning("Common-loading projection reached max_iter without convergence.", call. = FALSE)
  }
  if (alignment_method == "Procrustes") {
    failed <- which(!convergence$procrustes$converged)
    failed <- setdiff(failed, common$common_component$alignment$reference_group)
    if (length(failed)) {
      groups <- if (is.null(data_info$group_names)) failed else data_info$group_names[failed]
      warning("Procrustes alignment reached max_iter_procrustes without convergence for group(s): ",
              paste(groups, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (refit_method == "ALS" && !convergence$als$converged) {
    warning("ALS refitting reached max_iter_als without convergence.", call. = FALSE)
  }
  structure(list(call = call, dimensions = data_info, rank = rank, main_effect = main$main_effect,
                 common_component = common$common_component, check_Y = main$check_Y, convergence = convergence),
            class = "mmefm_fit")
}

#' Reconstruct fitted MMEFM components
#'
#' Recover the estimated signal or one of its global and local components from an MMEFM fit.
#' @param object An `mmefm_fit` returned by [est_MMEFM()].
#' @param component Component to reconstruct. `"all"` (default) sums all four structured components; `"main"` sums global and local main effects. `"global_main"` and `"local_main"` include their respective grand means, row main effects, and column main effects. `"global_common"` and `"local_common"` give the corresponding matrix-factor interactions.
#' @param ... Further arguments, currently unused.
#' @details Global main effects use `A1_tilde` and `B1_tilde` with the global `alpha_hat` and `beta_hat` scores. Local main effects use `A2_hat` and `B2_hat` with group-local scores. These final factor representations differ from the direct moment estimates in `object$main_effect$direct`.
#'
#' Global common components use `Q_tilde`, `G_tilde`, and `J_tilde`. Local common components use `R_hat`, `F_hat`, and `C_hat`. `check_Y` is an intermediate residual and is not an additional fitted component.
#' @return A list of `T x p_m x q_m` arrays preserving the data's group names and dimensions.
#' @seealso [est_MMEFM()] for an executable example, [residuals.mmefm_fit()], [summary.mmefm_fit()]
#' @export
fitted.mmefm_fit <- function(object, component = c("all", "main", "global_main", "local_main",
  "global_common", "local_common"), ...) {
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

#' Compute residuals of an MMEFM fit
#'
#' Subtract the complete final fitted signal from the observed group arrays.
#' @inheritParams fitted.mmefm_fit
#' @details The final fitted signal contains low-rank global and local main effects and common components. It differs from the direct additive estimate removed when constructing `object$check_Y`. Hence `check_Y` and the final model residual are distinct.
#'
#' The observations are recovered from the retained direct moments and `check_Y`; a separate observed-data copy is not needed.
#' @return A list of `T x p_m x q_m` residual arrays preserving group names and dimensions.
#' @seealso [fitted.mmefm_fit()], [est_MMEFM()] for an executable example
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

#' Summarize and print an MMEFM fit
#'
#' Report the fitted model's dimensions, global and local ranks, alignment method, and refitting method.
#' @inheritParams fitted.mmefm_fit
#' @details The summary omits rank-selection diagnostics and convergence information. Inspect `object$rank$diagnostics` and `object$convergence` for these, respectively. Printing a fit displays the same model summary.
#' @return `summary()` returns an object of class `summary.mmefm_fit` containing `call`, `dimensions`, the eight statistical `rank` values, `alignment_method`, and `refit_method`. The print methods return their input invisibly.
#' @seealso [est_MMEFM()] for an executable example, [fitted.mmefm_fit()], [residuals.mmefm_fit()]
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
