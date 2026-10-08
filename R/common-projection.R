# Directions are unit column vectors; check_Y has already passed the data contract.
.common_global_psd <- function(check_Y, directions, side) {
  side <- match.arg(side, c("row", "column"))
  M <- length(check_Y)
  TT <- dim(check_Y[[1]])[1L]
  projections <- vector("list", M)
  for (m in seq_len(M)) {
    p_m <- dim(check_Y[[m]])[2L]
    q_m <- dim(check_Y[[m]])[3L]
    if (side == "row") {
      projected <- matrix(check_Y[[m]], TT * p_m, q_m) %*% directions[[m]]
      projections[[m]] <- matrix(projected, TT, p_m)
    } else {
      projected <- matrix(aperm(check_Y[[m]], c(1L, 3L, 2L)), TT * q_m, p_m) %*% directions[[m]]
      projections[[m]] <- matrix(projected, TT, q_m)
    }
  }
  X <- vector("list", M)
  names(X) <- names(check_Y)
  for (m in seq_len(M)) {
    X[[m]] <- matrix(0, ncol(projections[[m]]), ncol(projections[[m]]))
    for (n in setdiff(seq_len(M), m)) {
      cross_proj <- crossprod(projections[[m]], projections[[n]])
      X[[m]] <- X[[m]] + tcrossprod(cross_proj)
    }
    X[[m]] <- (X[[m]] + t(X[[m]])) / (2 * TT^2 * (M - 1L))
  }
  X
}

# Cache candidate projections, then update groups sequentially from candidate one.
.select_common_directions <- function(check_Y, pools, side) {
  side <- match.arg(side, c("row", "column"))
  M <- length(check_Y)
  TT <- dim(check_Y[[1]])[1L]
  projections <- vector("list", M)
  selected_proj <- directions <- vector("list", M)
  names(directions) <- names(check_Y)
  for (m in seq_len(M)) {
    p_m <- dim(check_Y[[m]])[2L]
    q_m <- dim(check_Y[[m]])[3L]
    V <- do.call(cbind, pools[[m]])
    if (side == "row") {
      projected <- matrix(check_Y[[m]], TT * p_m, q_m) %*% V
      d <- p_m
    } else {
      projected <- matrix(aperm(check_Y[[m]], c(1L, 3L, 2L)), TT * q_m, p_m) %*% V
      d <- q_m
    }
    projections[[m]] <- vector("list", length(pools[[m]]))
    for (j in seq_along(pools[[m]])) projections[[m]][[j]] <- matrix(projected[, j, drop = FALSE], TT, d)
    selected_proj[[m]] <- projections[[m]][[1L]]
    directions[[m]] <- pools[[m]][[1L]]
  }
  for (m in seq_len(M)) {
    scores <- numeric(length(pools[[m]]))
    for (j in seq_along(pools[[m]])) {
      X <- matrix(0, ncol(projections[[m]][[j]]), ncol(projections[[m]][[j]]))
      for (n in setdiff(seq_len(M), m)) {
        cross_proj <- crossprod(projections[[m]][[j]], selected_proj[[n]])
        X <- X + tcrossprod(cross_proj)
      }
      # The common T^2(M-1) divisor does not affect candidate ordering.
      scores[j] <- eigen((X + t(X)) / 2, symmetric = TRUE, only.values = TRUE)$values[1L]
    }
    selected <- which.max(scores)
    selected_proj[[m]] <- projections[[m]][[selected]]
    directions[[m]] <- pools[[m]][[selected]]
  }
  directions
}

# Algorithm 1 updates automatic ranks inside the projection iteration; the caller owns the RNG seed.
.initial_common_loadings <- function(check_Y, kr = NULL, kc = NULL, K0 = 20L, max_iter = 20L, tol = 1e-4, verbose = FALSE) {
  if (is.null(kr) != is.null(kc)) stop("kr and kc must either both be supplied or both be NULL.")
  automatic_rank <- is.null(kr)
  M <- length(check_Y)
  pools <- vector("list", M)
  for (m in seq_len(M)) {
    q_m <- dim(check_Y[[m]])[3L]
    pools[[m]] <- vector("list", K0)
    for (k in seq_len(K0)) {
      v <- matrix(0, q_m, 1L)
      v[sample(q_m, floor(q_m / 2)), 1L] <- 1
      pools[[m]][[k]] <- v / sqrt(sum(v^2))
    }
  }
  directions <- .select_common_directions(check_Y, pools, "row")
  row_psd <- .common_global_psd(check_Y, directions, "row")
  if (automatic_rank) {
    row_rank_fit <- .paper_eigen_ratio(row_psd)
    kr <- row_rank_fit$rank
    Q_hat <- lapply(row_rank_fit$eigenvectors, function(V) V[, seq_len(kr), drop = FALSE])
  } else {
    Q_hat <- lapply(row_psd, .leading_eigenvectors, k = kr)
  }
  for (m in seq_len(M)) {
    pools[[m]] <- vector("list", kr)
    for (k in seq_len(kr)) pools[[m]][[k]] <- Q_hat[[m]][, k, drop = FALSE]
  }
  directions <- .select_common_directions(check_Y, pools, "column")
  column_psd <- .common_global_psd(check_Y, directions, "column")
  if (automatic_rank) {
    column_rank_fit <- .paper_eigen_ratio(column_psd)
    kc <- column_rank_fit$rank
    J_hat <- lapply(column_rank_fit$eigenvectors, function(V) V[, seq_len(kc), drop = FALSE])
  } else {
    J_hat <- lapply(column_psd, .leading_eigenvectors, k = kc)
  }

  iterations <- 0L
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    iterations <- iter
    Q_old <- Q_hat
    J_old <- J_hat
    previous_rank <- c(kr, kc)
    for (m in seq_len(M)) {
      pools[[m]] <- vector("list", kr)
      for (k in seq_len(kr)) pools[[m]][[k]] <- Q_hat[[m]][, k, drop = FALSE]
    }
    directions <- .select_common_directions(check_Y, pools, "column")
    column_psd <- .common_global_psd(check_Y, directions, "column")
    if (automatic_rank) {
      column_rank_fit <- .paper_eigen_ratio(column_psd)
      kc <- column_rank_fit$rank
      J_hat <- lapply(column_rank_fit$eigenvectors, function(V) V[, seq_len(kc), drop = FALSE])
    } else {
      J_hat <- lapply(column_psd, .leading_eigenvectors, k = kc)
    }
    for (m in seq_len(M)) {
      pools[[m]] <- vector("list", kc)
      for (k in seq_len(kc)) pools[[m]][[k]] <- J_hat[[m]][, k, drop = FALSE]
    }
    directions <- .select_common_directions(check_Y, pools, "row")
    row_psd <- .common_global_psd(check_Y, directions, "row")
    if (automatic_rank) {
      row_rank_fit <- .paper_eigen_ratio(row_psd)
      kr <- row_rank_fit$rank
      Q_hat <- lapply(row_rank_fit$eigenvectors, function(V) V[, seq_len(kr), drop = FALSE])
    } else {
      Q_hat <- lapply(row_psd, .leading_eigenvectors, k = kr)
    }
    discrepancy <- Inf
    if (identical(c(kr, kc), previous_rank)) {
      discrepancy <- 0
      for (m in seq_len(M)) discrepancy <- discrepancy + .subspace_distance(Q_hat[[m]], Q_old[[m]]) + .subspace_distance(J_hat[[m]], J_old[[m]])
    }
    if (verbose) cat(sprintf("projection iteration %d: discrepancy=%.5f\n", iter, discrepancy))
    if (is.finite(discrepancy) && discrepancy < tol) {
      converged <- TRUE
      break
    }
  }
  result <- list(Q_hat = Q_hat, J_hat = J_hat, iterations = iterations, converged = converged)
  if (automatic_rank) {
    result$kr <- kr
    result$kc <- kc
    result$row_rank_fit <- row_rank_fit
    result$column_rank_fit <- column_rank_fit
  }
  result
}
