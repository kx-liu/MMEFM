.local_common_psd <- function(check_Y_m, Q_hat_m, J_hat_m) {
  TT <- dim(check_Y_m)[1L]
  p <- dim(check_Y_m)[2L]
  q <- dim(check_Y_m)[3L]
  X_R <- matrix(0, p, p)
  X_C <- matrix(0, q, q)
  for (t in seq_len(TT)) {
    Yt <- matrix(check_Y_m[t, , , drop = FALSE], p, q)
    Y_perp_J <- Yt - Yt %*% J_hat_m %*% t(J_hat_m)
    Y_perp_Q <- Yt - Q_hat_m %*% crossprod(Q_hat_m, Yt)
    X_R <- X_R + tcrossprod(Y_perp_J)
    X_C <- X_C + crossprod(Y_perp_Q)
  }
  list(row = (X_R + t(X_R)) / (2 * TT), column = (X_C + t(X_C)) / (2 * TT))
}

# Shared least-squares G for arbitrary, possibly rescaled, global loadings.
.estimate_common_factors <- function(Y, Q, J) {
  TT <- dim(Y[[1L]])[1L]
  kr <- ncol(Q[[1L]])
  kc <- ncol(J[[1L]])
  gram <- matrix(0, kr * kc, kr * kc)
  numerator <- matrix(0, kr * kc, TT)
  for (m in seq_along(Y)) {
    gram <- gram + kronecker(crossprod(J[[m]]), crossprod(Q[[m]]))
    numerator <- numerator + .core_score_matrix(.project_common_core(Y[[m]], Q[[m]], J[[m]]))
  }
  scores <- .solve_gram(gram, numerator)
  aperm(array(scores, c(kr, kc, TT)), c(3L, 1L, 2L))
}

.refit_common_als <- function(Y, Q, J, G, max_iter_als, tol_als, verbose = FALSE) {
  TT <- dim(G)[1L]
  kr <- dim(G)[2L]
  kc <- dim(G)[3L]
  previous_objective <- 0
  for (m in seq_along(Y)) {
    for (t in seq_len(TT)) {
      Yt <- matrix(Y[[m]][t, , , drop = FALSE], nrow(Q[[m]]), nrow(J[[m]]))
      Gt <- matrix(G[t, , , drop = FALSE], kr, kc)
      previous_objective <- previous_objective + sum((Yt - Q[[m]] %*% Gt %*% t(J[[m]]))^2)
    }
  }
  iterations <- 0L
  converged <- FALSE
  for (iter in seq_len(max_iter_als)) {
    for (m in seq_along(Y)) {
      numerator_Q <- matrix(0, nrow(Q[[m]]), kr)
      gram_Q <- matrix(0, kr, kr)
      JJ <- crossprod(J[[m]])
      for (t in seq_len(TT)) {
        Yt <- matrix(Y[[m]][t, , , drop = FALSE], nrow(Q[[m]]), nrow(J[[m]]))
        Gt <- matrix(G[t, , , drop = FALSE], kr, kc)
        numerator_Q <- numerator_Q + Yt %*% J[[m]] %*% t(Gt)
        gram_Q <- gram_Q + Gt %*% JJ %*% t(Gt)
      }
      Q[[m]] <- t(.solve_gram(t(gram_Q), t(numerator_Q)))
      numerator_J <- matrix(0, nrow(J[[m]]), kc)
      gram_J <- matrix(0, kc, kc)
      QQ <- crossprod(Q[[m]])
      for (t in seq_len(TT)) {
        Yt <- matrix(Y[[m]][t, , , drop = FALSE], nrow(Q[[m]]), nrow(J[[m]]))
        Gt <- matrix(G[t, , , drop = FALSE], kr, kc)
        numerator_J <- numerator_J + t(Yt) %*% Q[[m]] %*% Gt
        gram_J <- gram_J + t(Gt) %*% QQ %*% Gt
      }
      J[[m]] <- t(.solve_gram(t(gram_J), t(numerator_J)))
    }
    # The deliberate public convention updates the shared G after each full Q/J sweep.
    G <- .estimate_common_factors(Y, Q, J)
    objective <- 0
    for (m in seq_along(Y)) {
      for (t in seq_len(TT)) {
        Yt <- matrix(Y[[m]][t, , , drop = FALSE], nrow(Q[[m]]), nrow(J[[m]]))
        Gt <- matrix(G[t, , , drop = FALSE], kr, kc)
        objective <- objective + sum((Yt - Q[[m]] %*% Gt %*% t(J[[m]]))^2)
      }
    }
    iterations <- iter
    relative_change <- abs(previous_objective - objective) / max(previous_objective, 1e-15)
    if (verbose) {
      cat(sprintf("ALS iteration %d: objective=%.6e, relative change=%.2e\n", iter, objective, relative_change))
    }
    if (is.finite(relative_change) && relative_change < tol_als) {
      converged <- TRUE
      break
    }
    previous_objective <- objective
  }
  list(Q_tilde = Q, J_tilde = J, G_tilde = G, iterations = iterations, converged = converged)
}

.estimate_common_component <- function(
    check_Y, rank, K0 = 20L, max_iter = 20L, tol = 1e-4,
    alignment_method = c("VanLoan", "Procrustes"), max_iter_procrustes = 20L, tol_procrustes = 1e-4,
    refit_method = c("VanLoan", "ALS"), max_iter_als = 20L, tol_als = 1e-4, verbose = FALSE) {
  alignment_method <- match.arg(alignment_method)
  refit_method <- match.arg(refit_method)
  initial <- .initial_common_loadings(check_Y, rank$kr, rank$kc, K0, max_iter, tol, verbose)
  M <- length(check_Y)
  TT <- dim(check_Y[[1L]])[1L]
  R_hat <- C_hat <- F_hat <- Y_no_local <- vector("list", M)
  names(R_hat) <- names(C_hat) <- names(F_hat) <- names(Y_no_local) <- names(check_Y)
  for (m in seq_len(M)) {
    p <- dim(check_Y[[m]])[2L]
    q <- dim(check_Y[[m]])[3L]
    Q <- initial$Q_hat[[m]]
    psd <- .local_common_psd(check_Y[[m]], Q, initial$J_hat[[m]])
    R_hat[[m]] <- .leading_eigenvectors(psd$row, rank$kr_m[m])
    C_hat[[m]] <- .leading_eigenvectors(psd$column, rank$kc_m[m])
    projected_R <- R_hat[[m]] - Q %*% crossprod(Q, R_hat[[m]])
    left_F <- .solve_gram(crossprod(R_hat[[m]], projected_R), t(projected_R))
    F_hat[[m]] <- array(0, c(TT, rank$kr_m[m], rank$kc_m[m]))
    Y_no_local[[m]] <- array(0, c(TT, p, q))
    for (t in seq_len(TT)) {
      Yt <- matrix(check_Y[[m]][t, , , drop = FALSE], p, q)
      Ft <- left_F %*% Yt %*% C_hat[[m]]
      F_hat[[m]][t, , ] <- Ft
      Y_no_local[[m]][t, , ] <- Yt - R_hat[[m]] %*% Ft %*% t(C_hat[[m]])
    }
  }
  aligned <- .align_common_coordinates(Y_no_local, initial$Q_hat, initial$J_hat,
                                       alignment_method, max_iter_procrustes, tol_procrustes)
  G_hat <- aligned$G_hat
  als <- NULL
  if (refit_method == "VanLoan") {
    Z_G <- .core_score_matrix(G_hat)
    gram <- tcrossprod(Z_G)
    Q_tilde <- J_tilde <- J_kron_Q_hat <- vector("list", M)
    names(Q_tilde) <- names(J_tilde) <- names(J_kron_Q_hat) <- names(check_Y)
    for (m in seq_len(M)) {
      p <- dim(check_Y[[m]])[2L]
      q <- dim(check_Y[[m]])[3L]
      response <- matrix(aperm(Y_no_local[[m]], c(2L, 3L, 1L)), p * q, TT)
      J_kron_Q_hat[[m]] <- t(.solve_gram(t(gram), t(tcrossprod(response, Z_G))))
      factors <- .van_loan_factors(J_kron_Q_hat[[m]], rank$kr, rank$kc, p, q)
      Q_tilde[[m]] <- factors$Q
      J_tilde[[m]] <- factors$J
    }
    aligned$alignment$J_kron_Q_hat <- J_kron_Q_hat
    G_tilde <- G_hat
  } else {
    refit <- .refit_common_als(Y_no_local, aligned$Q_aligned, aligned$J_aligned, G_hat,
                             max_iter_als, tol_als, verbose)
    Q_tilde <- refit$Q_tilde
    J_tilde <- refit$J_tilde
    G_tilde <- refit$G_tilde
    als <- list(iterations = refit$iterations, converged = refit$converged)
  }
  list(common_component = list(
         global = list(Q_hat = initial$Q_hat, J_hat = initial$J_hat, G_hat = G_hat,
                       Q_tilde = Q_tilde, J_tilde = J_tilde, G_tilde = G_tilde),
         local = list(R_hat = R_hat, C_hat = C_hat, F_hat = F_hat), alignment = aligned$alignment),
       convergence = list(projection = list(iterations = initial$iterations, converged = initial$converged),
                          procrustes = aligned$procrustes, als = als))
}
