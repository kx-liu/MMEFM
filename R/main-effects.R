# Step 1 retains time as the first axis, including when T = 1.
.main_effect_moments <- function(Xt) {
  M <- length(Xt)
  TT <- dim(Xt[[1]])[1L]
  c_hat <- a_hat <- b_hat <- check_Y <- vector("list", M)
  names(c_hat) <- names(a_hat) <- names(b_hat) <- names(check_Y) <- names(Xt)
  for (m in seq_len(M)) {
    p_m <- dim(Xt[[m]])[2L]
    q_m <- dim(Xt[[m]])[3L]
    c_hat[[m]] <- matrix(rowMeans(matrix(Xt[[m]], TT, p_m * q_m)), TT, 1L)
    a_hat[[m]] <- sweep(rowSums(Xt[[m]], dims = 2L) / q_m, 1L, c_hat[[m]][, 1L], "-")
    b_hat[[m]] <- sweep(rowSums(aperm(Xt[[m]], c(1L, 3L, 2L)), dims = 2L) / p_m, 1L, c_hat[[m]][, 1L], "-")
    check_Y[[m]] <- array(0, c(TT, p_m, q_m))
    for (t in seq_len(TT)) {
      Yt <- matrix(Xt[[m]][t, , , drop = FALSE], p_m, q_m)
      residual <- Yt - matrix(c_hat[[m]][t, 1L], p_m, q_m)
      residual <- sweep(residual, 1L, a_hat[[m]][t, ], "-")
      check_Y[[m]][t, , ] <- sweep(residual, 2L, b_hat[[m]][t, ], "-")
    }
  }
  list(direct = list(c_hat = c_hat, a_hat = a_hat, b_hat = b_hat), check_Y = check_Y)
}

# Cross-group covariance isolates the global main-effect loading spaces.
.main_effect_global_psd <- function(a_hat, b_hat, TT) {
  M <- length(a_hat)
  Omega_A1 <- Omega_B1 <- vector("list", M)
  names(Omega_A1) <- names(Omega_B1) <- names(a_hat)
  for (m in seq_len(M)) {
    Omega_A1[[m]] <- matrix(0, ncol(a_hat[[m]]), ncol(a_hat[[m]]))
    Omega_B1[[m]] <- matrix(0, ncol(b_hat[[m]]), ncol(b_hat[[m]]))
    for (n in setdiff(seq_len(M), m)) {
      S_A1_mn <- crossprod(a_hat[[m]], a_hat[[n]]) / TT
      S_B1_mn <- crossprod(b_hat[[m]], b_hat[[n]]) / TT
      Omega_A1[[m]] <- Omega_A1[[m]] + tcrossprod(S_A1_mn)
      Omega_B1[[m]] <- Omega_B1[[m]] + tcrossprod(S_B1_mn)
    }
    Omega_A1[[m]] <- Omega_A1[[m]] / (M - 1L)
    Omega_B1[[m]] <- Omega_B1[[m]] / (M - 1L)
  }
  list(row = Omega_A1, column = Omega_B1)
}

.main_effect_local_psd <- function(a_hat, b_hat, A1_hat, B1_hat, TT) {
  M <- length(a_hat)
  Omega_A2 <- Omega_B2 <- vector("list", M)
  names(Omega_A2) <- names(Omega_B2) <- names(a_hat)
  for (m in seq_len(M)) {
    P_A <- diag(ncol(a_hat[[m]])) - tcrossprod(A1_hat[[m]])
    P_B <- diag(ncol(b_hat[[m]])) - tcrossprod(B1_hat[[m]])
    S_A <- crossprod(a_hat[[m]]) / TT
    S_B <- crossprod(b_hat[[m]]) / TT
    Omega_A2[[m]] <- tcrossprod(S_A %*% P_A)
    Omega_B2[[m]] <- tcrossprod(S_B %*% P_B)
  }
  list(row = Omega_A2, column = Omega_B2)
}

# Xt and rank have passed the input contracts; lambda is the optional ridge extension.
.estimate_main_effect <- function(Xt, rank, lambda = 0) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite nonnegative number.")
  }
  M <- length(Xt)
  TT <- dim(Xt[[1]])[1L]
  moments <- .main_effect_moments(Xt)
  c_hat <- moments$direct$c_hat
  a_hat <- moments$direct$a_hat
  b_hat <- moments$direct$b_hat
  mu_hat <- Reduce("+", c_hat) / M
  local_mu_hat <- lapply(c_hat, function(cm) cm - mu_hat)

  A1_hat <- B1_hat <- A2_hat <- B2_hat <- vector("list", M)
  names(A1_hat) <- names(B1_hat) <- names(A2_hat) <- names(B2_hat) <- names(Xt)
  global_psd <- .main_effect_global_psd(a_hat, b_hat, TT)
  for (m in seq_len(M)) {
    A1_hat[[m]] <- .leading_eigenvectors(global_psd$row[[m]], rank$r1)
    B1_hat[[m]] <- .leading_eigenvectors(global_psd$column[[m]], rank$l1)
  }
  local_psd <- .main_effect_local_psd(a_hat, b_hat, A1_hat, B1_hat, TT)
  for (m in seq_len(M)) {
    A2_hat[[m]] <- .leading_eigenvectors(local_psd$row[[m]], rank$r2[m])
    B2_hat[[m]] <- .leading_eigenvectors(local_psd$column[[m]], rank$l2[m])
  }

  local_alpha_hat <- local_beta_hat <- y_hat <- w_hat <- vector("list", M)
  names(local_alpha_hat) <- names(local_beta_hat) <- names(y_hat) <- names(w_hat) <- names(Xt)
  alpha_projector <- beta_projector <- matrix(0, TT, TT)
  for (m in seq_len(M)) {
    P_A <- diag(ncol(a_hat[[m]])) - tcrossprod(A1_hat[[m]])
    P_B <- diag(ncol(b_hat[[m]])) - tcrossprod(B1_hat[[m]])
    P_A2 <- P_A %*% A2_hat[[m]]
    P_B2 <- P_B %*% B2_hat[[m]]
    local_alpha_hat[[m]] <- t(.solve_gram(crossprod(A2_hat[[m]], P_A2), crossprod(P_A2, t(a_hat[[m]])), lambda))
    local_beta_hat[[m]] <- t(.solve_gram(crossprod(B2_hat[[m]], P_B2), crossprod(P_B2, t(b_hat[[m]])), lambda))
    y_hat[[m]] <- a_hat[[m]] - local_alpha_hat[[m]] %*% t(A2_hat[[m]])
    w_hat[[m]] <- b_hat[[m]] - local_beta_hat[[m]] %*% t(B2_hat[[m]])
    Z_alpha <- crossprod(A1_hat[[m]], t(y_hat[[m]]))
    Z_beta <- crossprod(B1_hat[[m]], t(w_hat[[m]]))
    alpha_projector <- alpha_projector + crossprod(Z_alpha, .solve_gram(tcrossprod(Z_alpha), Z_alpha, lambda))
    beta_projector <- beta_projector + crossprod(Z_beta, .solve_gram(tcrossprod(Z_beta), Z_beta, lambda))
  }
  # The pooled temporal eigenvectors fix the paper's T-normalized global scores.
  alpha_hat <- sqrt(TT) * .leading_eigenvectors(alpha_projector / M, rank$r1)
  beta_hat <- sqrt(TT) * .leading_eigenvectors(beta_projector / M, rank$l1)
  A1_tilde <- B1_tilde <- vector("list", M)
  names(A1_tilde) <- names(B1_tilde) <- names(Xt)
  for (m in seq_len(M)) {
    if (lambda == 0) {
      A1_tilde[[m]] <- crossprod(y_hat[[m]], alpha_hat) / TT
      B1_tilde[[m]] <- crossprod(w_hat[[m]], beta_hat) / TT
    } else {
      A1_tilde[[m]] <- t(.solve_gram(crossprod(alpha_hat), crossprod(alpha_hat, y_hat[[m]]), lambda))
      B1_tilde[[m]] <- t(.solve_gram(crossprod(beta_hat), crossprod(beta_hat, w_hat[[m]]), lambda))
    }
  }
  list(main_effect = list(
    direct = moments$direct,
    global = list(mu_hat = mu_hat, A1_hat = A1_hat, alpha_hat = alpha_hat, A1_tilde = A1_tilde, B1_hat = B1_hat, beta_hat = beta_hat, B1_tilde = B1_tilde),
    local = list(mu_hat = local_mu_hat, A2_hat = A2_hat, alpha_hat = local_alpha_hat, B2_hat = B2_hat, beta_hat = local_beta_hat)
  ), check_Y = moments$check_Y)
}
