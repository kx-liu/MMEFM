# A single PSD supplies local ratios; a list supplies indexwise global maxima.
.paper_eigen_ratio <- function(S) {
  matrices <- if (is.matrix(S)) list(S) else S
  if (!is.list(matrices) || !length(matrices)) stop("Provide a PSD matrix or a nonempty list of PSD matrices.")
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
    if (any(values < -negative_roundoff)) stop("A rank-selection matrix is materially non-positive-semidefinite.")
    values[values < 0] <- 0
    eigenvalues[[m]] <- values
    eigenvectors[[m]] <- decomposition$vectors
  }
  d <- min(vapply(eigenvalues, length, integer(1L)))
  index <- seq_len(d - 1L)
  pooled_eigenvalues <- numeric(d)
  for (m in seq_along(matrices)) pooled_eigenvalues <- pmax(pooled_eigenvalues, eigenvalues[[m]][seq_len(d)])
  ratios <- rep(Inf, d - 1L)
  positive <- pooled_eigenvalues[index] > 0
  ratios[positive] <- pooled_eigenvalues[index[positive] + 1L] / pooled_eigenvalues[index[positive]]
  list(rank = as.integer(which.min(ratios)), eigenvalues = eigenvalues,
       pooled_eigenvalues = pooled_eigenvalues, ratios = ratios, eigenvectors = eigenvectors)
}

#' Select ranks for the Multilevel Main Effects Matrix Factor Model
#'
#' Estimates the eight main-effect and common-component loading ranks using the manuscript's unperturbed consecutive eigenvalue-ratio rule.
#'
#' @param Xt A list of at least two finite numeric arrays with dimensions `T x p_m x q_m`. All groups share `T >= 1`; row and column dimensions may differ between groups and must each be at least two. Group names are optional but must be nonempty and unique when supplied.
#' @param K0 Positive whole number of random half-panel column directions per group used to initialize Algorithm 1. Defaults to 20.
#' @param max_iter Positive whole number limiting Algorithm 1 update iterations. Defaults to 20.
#' @param tol Finite positive tolerance for the total loading-space discrepancy in Algorithm 1.
#' @param seed Nonnegative whole number within R's integer range controlling Algorithm 1's random directions. The same seed gives reproducible results under the same RNG kind; the caller's prior RNG state is restored on return or error.
#' @param verbose One nonmissing logical indicating whether to report Algorithm 1 progress.
#'
#' @details Global ranks use consecutive ratios of the indexwise maxima of group eigenvalue spectra, over the shared eigenvalue index range. Local ranks use the consecutive ratios of each group's own spectrum. All available consecutive ratios are considered, with a zero denominator assigned an infinite ratio. Ties select the first minimizing index. Ranks are positive in this version.
#'
#' Global common ranks `kr` and `kc` are re-estimated within Algorithm 1. Convergence requires unchanged ranks and a total loading-space discrepancy below `tol`; the final iterate is returned with a warning if `max_iter` is reached without convergence.
#'
#' Rank matrices are symmetrized before eigendecomposition. Only floating-point-scale negative eigenvalues are set to zero; materially indefinite matrices produce an error.
#'
#' @return An object of class `mmefm_rank` with scalar integer global ranks `r1`, `l1`, `kr`, and `kc`, and group-indexed integer vectors `r2`, `l2`, `kr_m`, and `kc_m`. Here `r1`/`r2` describe global/local row main effects, `l1`/`l2` global/local column main effects, and `kr`/`kr_m` and `kc`/`kc_m` global/local row and column common components. Group vectors retain the names and order of `Xt`.
#'
#' The additional `diagnostics` field has exactly the same eight rank names. Global entries contain `eigenvalues` (group spectra), `pooled_eigenvalues` (indexwise spectral maxima), and `ratios`. Local entries are group lists containing only ordinary `eigenvalues` and `ratios` vectors. Selected ranks, eigenvectors, and iteration traces are omitted.
#'
#' Under IC1, global/local row rank sums `r1 + r2[m]` and `kr + kr_m[m]` must be at most `p_m - 1`; column sums `l1 + l2[m]` and `kc + kc_m[m]` must be at most `q_m - 1`. Equality at these centered-subspace dimensions is allowed. Infeasible automatic selections error without capping or modifying ranks.
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
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be one finite positive number.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be one nonmissing logical.")

  moments <- .main_effect_moments(Xt)
  global_psd <- .main_effect_global_psd(moments$direct$a_hat, moments$direct$b_hat, data_info$T)
  r1_fit <- .paper_eigen_ratio(global_psd$row)
  l1_fit <- .paper_eigen_ratio(global_psd$column)
  A1_hat <- lapply(r1_fit$eigenvectors, function(V) V[, seq_len(r1_fit$rank), drop = FALSE])
  B1_hat <- lapply(l1_fit$eigenvectors, function(V) V[, seq_len(l1_fit$rank), drop = FALSE])
  local_psd <- .main_effect_local_psd(moments$direct$a_hat, moments$direct$b_hat, A1_hat, B1_hat, data_info$T)
  r2_fit <- l2_fit <- kr_m_fit <- kc_m_fit <- vector("list", data_info$M)
  names(r2_fit) <- names(l2_fit) <- names(kr_m_fit) <- names(kc_m_fit) <- data_info$group_names
  for (m in seq_len(data_info$M)) {
    r2_fit[[m]] <- .paper_eigen_ratio(local_psd$row[[m]])
    l2_fit[[m]] <- .paper_eigen_ratio(local_psd$column[[m]])
  }

  common <- .with_preserved_seed(seed,
    .initial_common_loadings(moments$check_Y, K0 = K0, max_iter = max_iter, tol = tol, verbose = verbose))
  if (!common$converged) warning("Rank-selection projection reached max_iter without convergence.", call. = FALSE)
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
