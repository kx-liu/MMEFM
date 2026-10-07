# Column-major vec ordering, retaining (kr * kc) x T even at rank one.
.core_score_matrix <- function(core) {
  matrix(aperm(core, c(2L, 3L, 1L)), prod(dim(core)[2:3]), dim(core)[1L])
}

.project_common_core <- function(Y, Q, J) {
  TT <- dim(Y)[1L]
  core <- array(0, c(TT, ncol(Q), ncol(J)))
  for (t in seq_len(TT)) {
    Yt <- matrix(Y[t, , , drop = FALSE], nrow(Q), nrow(J))
    core[t, , ] <- crossprod(Q, Yt) %*% J
  }
  core
}

# Rearrange the J-by-Q blocks so a Kronecker product becomes a rank-one matrix.
.van_loan_factors <- function(Xi, kr, kc, p = kr, q = kc) {
  rearranged <- matrix(0, p * kr, q * kc)
  for (j in seq_len(kc)) {
    for (i in seq_len(q)) {
      rows <- (i - 1L) * p + seq_len(p)
      cols <- (j - 1L) * kr + seq_len(kr)
      rearranged[, (j - 1L) * q + i] <- as.vector(Xi[rows, cols, drop = FALSE])
    }
  }
  triplet <- svd(rearranged, nu = 1L, nv = 1L)
  if (!is.finite(triplet$d[1L]) || triplet$d[1L] <= 0) {
    stop("Van Loan rearrangement has no finite positive leading singular value.")
  }
  list(Q = matrix(sqrt(triplet$d[1L]) * triplet$u[, 1L], p, kr),
       J = matrix(sqrt(triplet$d[1L]) * triplet$v[, 1L], q, kc))
}

.align_common_coordinates <- function(Y, Q_hat, J_hat, method, max_iter_procrustes, tol_procrustes) {
  M <- length(Y)
  TT <- dim(Y[[1L]])[1L]
  kr <- ncol(Q_hat[[1L]])
  kc <- ncol(J_hat[[1L]])
  cores <- vector("list", M)
  names(cores) <- names(Y)
  norms <- numeric(M)
  for (m in seq_len(M)) {
    cores[[m]] <- .project_common_core(Y[[m]], Q_hat[[m]], J_hat[[m]])
    norms[m] <- sum(cores[[m]]^2)
  }
  reference_group <- which.max(norms)
  if (!length(reference_group) || !is.finite(norms[reference_group]) || norms[reference_group] <= 0) {
    stop("Projected global cores have no finite positive reference norm.")
  }
  Q_aligned <- Q_hat
  J_aligned <- J_hat
  aligned_cores <- cores
  Xi_hat <- Xi_Q_hat <- Xi_J_hat <- O_Q_hat <- O_J_hat <- NULL
  procrustes <- NULL
  if (method == "VanLoan") {
    Xi_hat <- Xi_Q_hat <- Xi_J_hat <- vector("list", M)
    names(Xi_hat) <- names(Xi_Q_hat) <- names(Xi_J_hat) <- names(Y)
    Z_ref <- .core_score_matrix(cores[[reference_group]])
  } else {
    O_Q_hat <- O_J_hat <- vector("list", M)
    names(O_Q_hat) <- names(O_J_hat) <- names(Y)
    procrustes <- list(iterations = structure(integer(M), names = names(Y)), converged = structure(rep(FALSE, M), names = names(Y)))
    procrustes$converged[reference_group] <- TRUE
  }
  for (m in seq_len(M)) {
    if (method == "VanLoan") {
      if (m == reference_group) {
        Xi_hat[[m]] <- diag(kr * kc)
        factors <- list(Q = diag(kr), J = diag(kc))
      } else {
        Z_m <- .core_score_matrix(cores[[m]])
        Xi_hat[[m]] <- t(solve(t(tcrossprod(Z_m)), t(tcrossprod(Z_ref, Z_m))))
        factors <- .van_loan_factors(Xi_hat[[m]], kr, kc)
      }
      Xi_Q_hat[[m]] <- factors$Q
      Xi_J_hat[[m]] <- factors$J
      Q_aligned[[m]] <- t(solve(t(factors$Q), t(Q_hat[[m]])))
      J_aligned[[m]] <- t(solve(t(factors$J), t(J_hat[[m]])))
      for (t in seq_len(TT)) {
        core_t <- matrix(cores[[m]][t, , , drop = FALSE], kr, kc)
        aligned_cores[[m]][t, , ] <- factors$Q %*% core_t %*% t(factors$J)
      }
    } else {
      O_Q <- diag(kr)
      O_J <- diag(kc)
      previous_score <- -Inf
      if (m != reference_group) {
        for (iter in seq_len(max_iter_procrustes)) {
          Psi_Q <- matrix(0, kr, kr)
          for (t in seq_len(TT)) {
            target <- matrix(cores[[m]][t, , , drop = FALSE], kr, kc)
            reference <- matrix(cores[[reference_group]][t, , , drop = FALSE], kr, kc)
            Psi_Q <- Psi_Q + target %*% O_J %*% t(reference)
          }
          decomposition <- svd(Psi_Q)
          O_Q <- decomposition$u %*% t(decomposition$v)
          Psi_J <- matrix(0, kc, kc)
          for (t in seq_len(TT)) {
            target <- matrix(cores[[m]][t, , , drop = FALSE], kr, kc)
            reference <- matrix(cores[[reference_group]][t, , , drop = FALSE], kr, kc)
            Psi_J <- Psi_J + t(target) %*% O_Q %*% reference
          }
          decomposition <- svd(Psi_J)
          O_J <- decomposition$u %*% t(decomposition$v)
          score <- 0
          for (t in seq_len(TT)) {
            target <- matrix(cores[[m]][t, , , drop = FALSE], kr, kc)
            reference <- matrix(cores[[reference_group]][t, , , drop = FALSE], kr, kc)
            score <- score + sum(reference * (t(O_Q) %*% target %*% O_J))
          }
          procrustes$iterations[m] <- iter
          # This floor belongs only to the established relative convergence criterion, not a solve.
          relative_change <- abs(score - previous_score) / max(abs(previous_score), 1e-12)
          if (is.finite(relative_change) && relative_change < tol_procrustes) {
            procrustes$converged[m] <- TRUE
            break
          }
          previous_score <- score
        }
      }
      O_Q_hat[[m]] <- O_Q
      O_J_hat[[m]] <- O_J
      Q_aligned[[m]] <- Q_hat[[m]] %*% O_Q
      J_aligned[[m]] <- J_hat[[m]] %*% O_J
      for (t in seq_len(TT)) {
        core_t <- matrix(cores[[m]][t, , , drop = FALSE], kr, kc)
        aligned_cores[[m]][t, , ] <- t(O_Q) %*% core_t %*% O_J
      }
    }
  }
  list(Q_aligned = Q_aligned, J_aligned = J_aligned, G_hat = Reduce("+", aligned_cores) / M,
       alignment = list(method = method, reference_group = reference_group,
                        Xi_hat = Xi_hat, Xi_Q_hat = Xi_Q_hat, Xi_J_hat = Xi_J_hat,
                        O_Q_hat = O_Q_hat, O_J_hat = O_J_hat, J_kron_Q_hat = NULL),
       procrustes = procrustes)
}
