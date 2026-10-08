# A single PSD supplies local ratios; a list supplies indexwise global maxima.
.paper_eigen_ratio <- function(S) {
  matrices <- if (is.matrix(S)) list(S) else S
  if (!is.list(matrices) || !length(matrices)) {
    stop("Provide a PSD matrix or a nonempty list of PSD matrices.")
  }
  eigenvalues <- eigenvectors <- vector("list", length(matrices))
  names(eigenvalues) <- names(eigenvectors) <- names(matrices)
  for (m in seq_along(matrices)) {
    Sm <- matrices[[m]]
    if (!is.matrix(Sm) || !is.numeric(Sm) || nrow(Sm) != ncol(Sm) ||
        nrow(Sm) < 2L || any(!is.finite(Sm))) {
      stop("Each PSD matrix must be finite, numeric, square, and have dimension at least two.")
    }
    decomposition <- eigen((Sm + t(Sm)) / 2, symmetric = TRUE)
    values <- decomposition$values
    # Backward error scales with matrix dimension and the spectrum, including at zero scale.
    negative_roundoff <- nrow(Sm) * .Machine$double.eps * max(abs(values))
    if (any(values < -negative_roundoff)) {
      stop("A rank-selection matrix is materially non-positive-semidefinite.")
    }
    values[values < 0] <- 0
    eigenvalues[[m]] <- values
    eigenvectors[[m]] <- decomposition$vectors
  }
  d <- min(vapply(eigenvalues, length, integer(1L)))
  index <- seq_len(d - 1L)
  # Pool eigenvalues indexwise across groups before taking consecutive ratios.
  pooled_eigenvalues <- numeric(d)
  for (m in seq_along(matrices)) {
    pooled_eigenvalues <- pmax(pooled_eigenvalues, eigenvalues[[m]][seq_len(d)])
  }
  ratios <- rep(Inf, d - 1L)
  positive <- pooled_eigenvalues[index] > 0
  ratios[positive] <- pooled_eigenvalues[index[positive] + 1L] / pooled_eigenvalues[index[positive]]
  list(rank = as.integer(which.min(ratios)), eigenvalues = eigenvalues,
       pooled_eigenvalues = pooled_eigenvalues, ratios = ratios, eigenvectors = eigenvectors)
}

#' Select global and local MMEFM ranks
#'
#' Estimate the loading ranks for global and local row main effects, column main effects, and common components. Selection uses the manuscript's unperturbed consecutive eigenvalue-ratio rule.
#'
#' @param Xt A list of at least two finite numeric arrays with dimensions `T x p_m x q_m`, in time, row, and column order. All groups share `T >= 1`; each spatial dimension is at least two. Group names are optional, but must be nonempty and unique when supplied.
#' @param K0 Positive whole number of random half-panel column directions per group for initializing the iterative projection procedure (Algorithm 1).
#' @param max_iter Positive whole-number limit on Algorithm 1 update iterations.
#' @param tol Finite positive tolerance for the sum of loading-space discrepancies in Algorithm 1.
#' @param seed Nonnegative whole number within R's integer range controlling random projection directions. Results are reproducible under the same seed and RNG kind. The caller's RNG state is restored on return or error.
#' @param verbose One nonmissing logical indicating whether to print Algorithm 1 progress.
#' @details
#' The eight ranks correspond to global/local row main effects (`r1`, `r2`), column main effects (`l1`, `l2`), row common loadings (`kr`, `kr_m`), and column common loadings (`kc`, `kc_m`). Ranks are positive in this version.
#'
#' For a local rank, the selected index minimizes \eqn{\lambda_{j+1}/\lambda_j} for the group's relevant moment matrix. For a global rank, eigenvalues are first pooled by taking their maximum across groups separately at each index. Consecutive ratios of that pooled spectrum are then minimized. The search uses every available ratio: through the smallest group dimension minus one for global ranks, and the group's dimension minus one for local ranks. A zero denominator gives an infinite ratio; ties select the first minimizing index.
#'
#' Global main-effect ranks use cross-group moments. Local main-effect ranks use moments projected off the estimated global spaces. Global common ranks are updated within Algorithm 1; local common ranks use residual moments after projection off the global common spaces. Projection convergence requires unchanged global common ranks and a loading-space discrepancy below `tol`. If `max_iter` is reached, the final iterate is used with a warning.
#'
#' The ranks are selected independently and then checked jointly. Under IC1, `r1 + r2[m]` and `kr + kr_m[m]` must not exceed `p_m - 1`; `l1 + l2[m]` and `kc + kc_m[m]` must not exceed `q_m - 1`. Equality is allowed. Infeasible selections error without capping or modifying ranks. Small-sample selection is not guaranteed to recover the true ranks.
#'
#' Integer controls must be within R's integer range. Moment matrices are symmetrized before eigendecomposition. Roundoff-size negative eigenvalues are set to zero; materially indefinite matrices error. No eigenvalue perturbation is added.
#' @return An `mmefm_rank` object containing the eight ranks and `diagnostics`. Global ranks `r1`, `l1`, `kr`, and `kc` are scalar integers. Local ranks `r2`, `l2`, `kr_m`, and `kc_m` are integer vectors retaining group names and order.
#'
#' `diagnostics` has the same eight names. Global entries contain group `eigenvalues`, indexwise `pooled_eigenvalues`, and `ratios`. Local entries are group lists containing `eigenvalues` and `ratios`. Eigenvectors and iteration traces are not returned. The rank object can be supplied directly to [est_MMEFM()].
#' @seealso [est_MMEFM()], [gen_MMEFM()]
#' @examples
#' # Exact additive and common signals with orthogonal temporal scores.
#' H <- matrix(1, 1, 1)
#' for (i in seq_len(4L)) {
#'   H <- rbind(cbind(H, H), cbind(H, -H))
#' }
#' u <- c(1, 1, -2)
#' v <- c(1, -1, 0)
#' Xt <- setNames(vector('list', 3L), c('a','b','c'))
#' for (m in seq_along(Xt)) {
#'   Xt[[m]] <- array(0, c(16L, 3L, 3L))
#'   for (t in seq_len(16L)) {
#'     row <- 3 * (m + 1) * H[t, 2L] * u + H[t, 2L + m] * v
#'     column <- 2 * (m + 2) * H[t, 6L] * u + H[t, 6L + m] * v
#'     Xt[[m]][t, , ] <- outer(row, rep(1, 3)) + outer(rep(1, 3), column) +
#'       5 * (m + 1) * (m + 2) * H[t, 10L] * outer(u, u) + H[t, 10L + m] * outer(v, v)
#'   }
#' }
#' selected <- select_MMEFM_rank(Xt, K0 = 8L, max_iter = 50L, seed = 2026L)
#' selected[c("r1", "l1", "kr", "kc")]
#' selected$diagnostics$r1$ratios
#' @export
select_MMEFM_rank <- function(Xt, K0 = 20L, max_iter = 20L, tol = 1e-4, seed = 2026L, verbose = FALSE) {
  data_info <- .validate_Xt(Xt)
  controls <- list(K0 = K0, max_iter = max_iter, seed = seed)
  for (name in names(controls)) {
    value <- controls[[name]]
    minimum <- if (name == "seed") 0L else 1L
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum || value > .Machine$integer.max) {
      stop(sprintf("%s must be one finite whole number >= %d within R's integer range.", name, minimum))
    }
  }
  K0 <- as.integer(K0)
  max_iter <- as.integer(max_iter)
  seed <- as.integer(seed)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one nonmissing logical.")
  }

  # Cross-group main-effect moments identify the global loading spaces.
  moments <- .main_effect_moments(Xt)
  global_psd <- .main_effect_global_psd(moments$direct$a_hat, moments$direct$b_hat, data_info$T)
  r1_fit <- .paper_eigen_ratio(global_psd$row)
  l1_fit <- .paper_eigen_ratio(global_psd$column)
  A1_hat <- lapply(r1_fit$eigenvectors, function(V) V[, seq_len(r1_fit$rank), drop = FALSE])
  B1_hat <- lapply(l1_fit$eigenvectors, function(V) V[, seq_len(l1_fit$rank), drop = FALSE])

  # Remove the global spaces before selecting group-local main-effect ranks.
  local_psd <- .main_effect_local_psd(moments$direct$a_hat, moments$direct$b_hat, A1_hat, B1_hat, data_info$T)
  r2_fit <- l2_fit <- kr_m_fit <- kc_m_fit <- vector("list", data_info$M)
  names(r2_fit) <- names(l2_fit) <- names(kr_m_fit) <- names(kc_m_fit) <- data_info$group_names
  for (m in seq_len(data_info$M)) {
    r2_fit[[m]] <- .paper_eigen_ratio(local_psd$row[[m]])
    l2_fit[[m]] <- .paper_eigen_ratio(local_psd$column[[m]])
  }

  common <- .with_preserved_seed(seed,
    .initial_common_loadings(moments$check_Y, K0 = K0, max_iter = max_iter, tol = tol, verbose = verbose))
  if (!common$converged) {
    warning("Rank-selection projection reached max_iter without convergence.", call. = FALSE)
  }

  # Project out global common spaces for each group-local common rank.
  for (m in seq_len(data_info$M)) {
    local_common <- .local_common_psd(moments$check_Y[[m]], common$Q_hat[[m]], common$J_hat[[m]])
    kr_m_fit[[m]] <- .paper_eigen_ratio(local_common$row)
    kc_m_fit[[m]] <- .paper_eigen_ratio(local_common$column)
  }
  ranks <- list(r1 = r1_fit$rank, l1 = l1_fit$rank,
                r2 = vapply(r2_fit, function(fit) fit$rank, integer(1L)),
                l2 = vapply(l2_fit, function(fit) fit$rank, integer(1L)),
                kr = common$kr, kc = common$kc,
                kr_m = vapply(kr_m_fit, function(fit) fit$rank, integer(1L)),
                kc_m = vapply(kc_m_fit, function(fit) fit$rank, integer(1L)))
  # Independently selected ranks must jointly fit the IC1-centered spaces.
  ranks <- .validate_rank(ranks, data_info)
  diagnostics <- list(r1 = r1_fit, l1 = l1_fit, r2 = r2_fit, l2 = l2_fit,
                      kr = common$row_rank_fit, kc = common$column_rank_fit, kr_m = kr_m_fit, kc_m = kc_m_fit)
  for (name in c("r1", "l1", "kr", "kc")) {
    diagnostics[[name]] <- diagnostics[[name]][c("eigenvalues", "pooled_eigenvalues", "ratios")]
  }
  for (name in c("r2", "l2", "kr_m", "kc_m")) {
    for (m in seq_len(data_info$M)) {
      fit <- diagnostics[[name]][[m]]
      diagnostics[[name]][[m]] <- list(eigenvalues = fit$eigenvalues[[1L]], ratios = fit$ratios)
    }
  }
  .new_mmefm_rank(ranks, diagnostics)
}
