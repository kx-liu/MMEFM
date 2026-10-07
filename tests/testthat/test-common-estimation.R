# Local factors use mutually orthogonal temporal directions, eliminating cross-group local moments.
.common_estimation_fixture <- function(kr = 1L, kc = 1L) {
  TT <- 12L
  p <- c(5L, 6L, 7L)
  q <- c(6L, 7L, 5L)
  temporal <- unname(contr.helmert(TT))
  temporal <- sweep(temporal, 2L, sqrt(colSums(temporal^2)), "/")
  Q <- J <- R <- C <- Y <- local <- global <- vector("list", 3L)
  names(Q) <- names(J) <- names(R) <- names(C) <- names(Y) <- names(local) <- names(global) <- c("a", "b", "c")
  G <- array(0, c(TT, kr, kc))
  for (t in seq_len(TT)) G[t, , ] <- matrix(c(6, 4, 3, 2)[seq_len(kr * kc)] * temporal[t, seq_len(kr * kc)], kr, kc)
  for (m in seq_len(3L)) {
    H_p <- unname(contr.helmert(p[m]))
    H_q <- unname(contr.helmert(q[m]))
    H_p <- sweep(H_p, 2L, sqrt(colSums(H_p^2)), "/")
    H_q <- sweep(H_q, 2L, sqrt(colSums(H_q^2)), "/")
    Q[[m]] <- (m + 1) * H_p[, seq_len(kr), drop = FALSE]
    J[[m]] <- (m + 2) * H_q[, seq_len(kc), drop = FALSE]
    R[[m]] <- H_p[, kr + 1L, drop = FALSE]
    C[[m]] <- H_q[, kc + 1L, drop = FALSE]
    Y[[m]] <- local[[m]] <- global[[m]] <- array(0, c(TT, p[m], q[m]))
    for (t in seq_len(TT)) {
      global[[m]][t, , ] <- Q[[m]] %*% matrix(G[t, , , drop = FALSE], kr, kc) %*% t(J[[m]])
      local[[m]][t, , ] <- 2 * temporal[t, kr * kc + m] * tcrossprod(R[[m]], C[[m]])
    }
    Y[[m]] <- global[[m]] + local[[m]]
  }
  list(Y = Y, Q = Q, J = J, R = R, C = C, G = G, local = local, global = global,
       rank = list(kr = kr, kc = kc, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L)))
}

.common_test_product <- function(Q, G, J) {
  out <- array(0, c(dim(G)[1L], nrow(Q), nrow(J)))
  for (t in seq_len(dim(G)[1L])) {
    out[t, , ] <- Q %*% matrix(G[t, , , drop = FALSE], ncol(Q), ncol(J)) %*% t(J)
  }
  out
}

test_that("local PSDs are direct projected Grams, including large removable signals", {
  Q <- diag(3L)[, 1L, drop = FALSE]
  J <- diag(4L)[, 1L, drop = FALSE]
  for (scale in c(1, 1e12)) {
    Y <- array(sin(seq_len(24L)), c(2L, 3L, 4L))
    Y[, 1L, 1L] <- scale
    psd <- .local_common_psd(Y, Q, J)
    X_R <- matrix(0, 3L, 3L)
    X_C <- matrix(0, 4L, 4L)
    for (t in seq_len(2L)) {
      Yt <- matrix(Y[t, , , drop = FALSE], 3L, 4L)
      row_projected <- Yt %*% (diag(4L) - tcrossprod(J))
      column_projected <- (diag(3L) - tcrossprod(Q)) %*% Yt
      X_R <- X_R + row_projected %*% t(row_projected) / 2
      X_C <- X_C + t(column_projected) %*% column_projected / 2
    }
    expect_equal(psd$row, X_R, tolerance = 1e-12)
    expect_equal(psd$column, X_C, tolerance = 1e-12)
    for (S in psd) {
      expect_identical(S, t(S))
      expect_true(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > -1e-12)
    }
  }
})

test_that("the complete estimator preserves rank-one shapes and exact global/local products", {
  for (kr in c(1L, 2L)) {
    fixture <- .common_estimation_fixture(kr, kr)
    set.seed(2026L)
    initial <- .initial_common_loadings(fixture$Y, kr, kr, K0 = 8L)
    set.seed(2026L)
    default <- .estimate_common_component(fixture$Y, fixture$rank, K0 = 8L)
    expect_identical(default$common_component$alignment$method, "VanLoan")
    expect_null(default$convergence$procrustes)
    expect_null(default$convergence$als)
    for (alignment in c("VanLoan", "Procrustes")) {
      for (refit in c("VanLoan", "ALS")) {
        set.seed(2026L)
        result <- .estimate_common_component(fixture$Y, fixture$rank, K0 = 8L,
                                             alignment_method = alignment, refit_method = refit)
        fit <- result$common_component
        expect_identical(names(result), c("common_component", "convergence"))
        expect_identical(names(fit$global), c("Q_hat", "J_hat", "G_hat", "Q_tilde", "J_tilde", "G_tilde"))
        expect_identical(names(fit$local), c("R_hat", "C_hat", "F_hat"))
        expect_identical(fit$global$Q_hat, initial$Q_hat)
        expect_identical(fit$global$J_hat, initial$J_hat)
        expect_identical(dim(fit$global$G_hat), c(12L, kr, kr))
        expect_identical(dim(fit$global$G_tilde), c(12L, kr, kr))
        expect_identical(fit$alignment$reference_group, 3L)
        for (field in c("Q_hat", "J_hat", "Q_tilde", "J_tilde")) expect_identical(names(fit$global[[field]]), names(fixture$Y))
        for (field in names(fit$local)) expect_identical(names(fit$local[[field]]), names(fixture$Y))
        for (m in seq_along(fixture$Y)) {
          expect_equal(tcrossprod(fit$local$R_hat[[m]]), tcrossprod(fixture$R[[m]]), tolerance = 1e-10)
          expect_equal(tcrossprod(fit$local$C_hat[[m]]), tcrossprod(fixture$C[[m]]), tolerance = 1e-10)
          local_fit <- .common_test_product(fit$local$R_hat[[m]], fit$local$F_hat[[m]], fit$local$C_hat[[m]])
          expect_equal(local_fit, fixture$local[[m]], tolerance = 1e-10)
          expect_equal(fixture$Y[[m]] - local_fit, fixture$global[[m]], tolerance = 1e-10)
          global_fit <- .common_test_product(fit$global$Q_tilde[[m]], fit$global$G_tilde, fit$global$J_tilde[[m]])
          expect_equal(global_fit, fixture$global[[m]], tolerance = 1e-10)
          expect_identical(dim(fit$global$Q_tilde[[m]]), c(dim(fixture$Y[[m]])[2L], kr))
          expect_identical(dim(fit$global$J_tilde[[m]]), c(dim(fixture$Y[[m]])[3L], kr))
          expect_identical(dim(fit$local$F_hat[[m]]), c(12L, 1L, 1L))
        }
        fields <- if (alignment == "VanLoan") c("Xi_hat", "Xi_Q_hat", "Xi_J_hat") else c("O_Q_hat", "O_J_hat")
        empty <- if (alignment == "VanLoan") c("O_Q_hat", "O_J_hat") else c("Xi_hat", "Xi_Q_hat", "Xi_J_hat")
        for (field in fields) expect_identical(names(fit$alignment[[field]]), names(fixture$Y))
        for (field in empty) expect_null(fit$alignment[[field]])
        if (alignment == "VanLoan") expect_null(result$convergence$procrustes) else expect_true(all(result$convergence$procrustes$converged))
        if (refit == "VanLoan") {
          expect_identical(fit$global$G_tilde, fit$global$G_hat)
          expect_identical(names(fit$alignment$J_kron_Q_hat), names(fixture$Y))
          for (m in seq_along(fixture$Y)) {
            expect_identical(dim(fit$alignment$J_kron_Q_hat[[m]]), as.integer(c(prod(dim(fixture$Y[[m]])[2:3]), kr^2L)))
          }
          expect_null(result$convergence$als)
        } else {
          expect_null(fit$alignment$J_kron_Q_hat)
          expect_true(is.logical(result$convergence$als$converged))
        }
      }
    }
    # The only random draws are those of Algorithm 1.
    set.seed(2026L)
    invisible(.initial_common_loadings(fixture$Y, kr, kr, K0 = 8L))
    projection_rng <- .Random.seed
    set.seed(2026L)
    invisible(.estimate_common_component(fixture$Y, fixture$rank, K0 = 8L))
    expect_identical(.Random.seed, projection_rng)
  }
})

test_that("alignment recovers two-sided coordinates and rejects singular score Grams", {
  for (kr in c(1L, 2L)) {
    fixture <- .common_estimation_fixture(kr, kr)
    Q <- J <- fixture$Q
    Y <- fixture$global
    for (m in seq_along(Y)) {
      Q[[m]] <- fixture$Q[[m]] / (m + 1)
      J[[m]] <- fixture$J[[m]] / (m + 2)
    }
    for (method in c("VanLoan", "Procrustes")) {
      Y_test <- Y
      if (method == "Procrustes") {
        # Remove scale differences so the groups differ only in orthogonal coordinates.
        Y_test <- lapply(seq_along(Y), function(m) Y[[m]] / ((m + 1) * (m + 2)))
        names(Y_test) <- names(Y)
        rotation <- if (kr == 1L) matrix(-1, 1L, 1L) else matrix(c(0, 1, -1, 0), 2L, 2L)
        Q[[1L]] <- Q[[1L]] %*% rotation
        J[[2L]] <- J[[2L]] %*% rotation
      }
      aligned <- .align_common_coordinates(Y_test, Q, J, method, 20L, 1e-10)
      ref <- aligned$alignment$reference_group
      reference <- .project_common_core(Y_test[[ref]], Q[[ref]], J[[ref]])
      for (m in seq_along(Y_test)) {
        core <- .project_common_core(Y_test[[m]], Q[[m]], J[[m]])
        transformed <- core
        if (method == "Procrustes") {
          O_Q <- aligned$alignment$O_Q_hat[[m]]
          O_J <- aligned$alignment$O_J_hat[[m]]
          expect_equal(crossprod(O_Q), diag(kr), tolerance = 1e-12)
          expect_equal(crossprod(O_J), diag(kr), tolerance = 1e-12)
        }
        for (t in seq_len(dim(core)[1L])) {
          core_t <- matrix(core[t, , , drop = FALSE], kr, kr)
          if (method == "VanLoan") {
            transformed[t, , ] <- aligned$alignment$Xi_Q_hat[[m]] %*% core_t %*% t(aligned$alignment$Xi_J_hat[[m]])
          } else {
            transformed[t, , ] <- t(O_Q) %*% core_t %*% O_J
          }
        }
        expect_equal(transformed, reference, tolerance = 1e-10)
        expect_equal(.common_test_product(aligned$Q_aligned[[m]], transformed, aligned$J_aligned[[m]]), Y_test[[m]], tolerance = 1e-10)
      }
      if (method == "VanLoan") expect_identical(ref, 3L)
      if (method == "Procrustes") {
        repeated <- .align_common_coordinates(Y_test, Q, J, method, 20L, 1e-10)
        expect_identical(aligned, repeated)
      }
    }
  }
  Y <- list(a = array(0, c(3L, 2L, 2L)), b = array(0, c(3L, 2L, 2L)))
  expect_error(.align_common_coordinates(Y, list(diag(2), diag(2)), list(diag(2), diag(2)), "VanLoan", 20L, 1e-4), "reference norm")
  Y$a[, 1L, 1L] <- 2
  Y$b[, 1L, 1L] <- 1
  expect_error(.align_common_coordinates(Y, list(diag(2), diag(2)), list(diag(2), diag(2)), "VanLoan", 20L, 1e-4), "singular")
  expect_identical(dim(.core_score_matrix(array(1:3, c(3L, 1L, 1L)))), c(1L, 3L))
})

test_that("ALS uses a shared least-squares G update after a full loading sweep", {
  fixture <- .common_estimation_fixture(2L, 2L)
  Y <- fixture$global
  # Deterministic perturbations make the initial shared fit inexact and the G update identifiable.
  Q <- fixture$Q
  J <- fixture$J
  for (m in seq_along(Y)) {
    Q[[m]] <- Q[[m]] + 0.07 * matrix(sin(seq_len(length(Q[[m]]))), nrow(Q[[m]]), 2L)
    J[[m]] <- J[[m]] + 0.05 * matrix(cos(seq_len(length(J[[m]]))), nrow(J[[m]]), 2L)
    Y[[m]] <- Y[[m]] + 0.02 * array(sin(seq_len(length(Y[[m]]))), dim(Y[[m]]))
  }
  G <- 0.9 * fixture$G
  initial_objective <- sum(vapply(seq_along(Y), function(m) sum((Y[[m]] - .common_test_product(Q[[m]], G, J[[m]]))^2), numeric(1L)))
  Q_expected <- Q
  J_expected <- J
  for (m in seq_along(Y)) {
    NumQ <- matrix(0, nrow(Q[[m]]), 2L)
    DenQ <- matrix(0, 2L, 2L)
    for (t in seq_len(12L)) {
      Gt <- matrix(G[t, , ], 2L, 2L)
      Yt <- matrix(Y[[m]][t, , ], nrow(Q[[m]]), nrow(J[[m]]))
      NumQ <- NumQ + Yt %*% J[[m]] %*% t(Gt)
      DenQ <- DenQ + Gt %*% crossprod(J[[m]]) %*% t(Gt)
    }
    Q_expected[[m]] <- NumQ %*% solve(DenQ)
    NumJ <- matrix(0, nrow(J[[m]]), 2L)
    DenJ <- matrix(0, 2L, 2L)
    for (t in seq_len(12L)) {
      Gt <- matrix(G[t, , ], 2L, 2L)
      Yt <- matrix(Y[[m]][t, , ], nrow(Q[[m]]), nrow(J[[m]]))
      NumJ <- NumJ + t(Yt) %*% Q_expected[[m]] %*% Gt
      DenJ <- DenJ + t(Gt) %*% crossprod(Q_expected[[m]]) %*% Gt
    }
    J_expected[[m]] <- NumJ %*% solve(DenJ)
  }
  design <- do.call(rbind, lapply(seq_along(Y), function(m) kronecker(J_expected[[m]], Q_expected[[m]])))
  expected_G <- array(0, c(12L, 2L, 2L))
  for (t in seq_len(12L)) {
    response <- unlist(lapply(Y, function(Ym) as.vector(Ym[t, , ])))
    expected_G[t, , ] <- solve(crossprod(design), crossprod(design, response))
  }
  fit <- .refit_common_als(Y, Q, J, G, 1L, 1e-4)
  expect_equal(fit$Q_tilde, Q_expected, tolerance = 1e-12)
  expect_equal(fit$J_tilde, J_expected, tolerance = 1e-12)
  expect_equal(fit$G_tilde, expected_G, tolerance = 1e-12)
  expect_gt(max(abs(fit$G_tilde - G)), 1e-6)
  expect_identical(fit$iterations, 1L)
  final_objective <- sum(vapply(seq_along(Y), function(m) sum((Y[[m]] - .common_test_product(fit$Q_tilde[[m]], fit$G_tilde, fit$J_tilde[[m]]))^2), numeric(1L)))
  expect_true(is.finite(final_objective))
  expect_lte(final_objective, initial_objective + 1e-12)
  completed <- .refit_common_als(Y, Q, J, G, 20L, 1e-4)
  completed_objective <- sum(vapply(seq_along(Y), function(m) sum((Y[[m]] - .common_test_product(completed$Q_tilde[[m]], completed$G_tilde, completed$J_tilde[[m]]))^2), numeric(1L)))
  expect_lte(completed_objective, final_objective + 1e-12)
})


test_that("rectangular global ranks retain vectorization and factor dimensions", {
  for (ranks in list(c(2L, 1L), c(1L, 2L))) {
    kr <- ranks[1L]
    kc <- ranks[2L]
    fixture <- .common_estimation_fixture(kr, kc)
    for (alignment in c("VanLoan", "Procrustes")) {
      for (refit in c("VanLoan", "ALS")) {
        set.seed(2026L)
        result <- .estimate_common_component(fixture$Y, fixture$rank, K0 = 8L,
                                             alignment_method = alignment, refit_method = refit)
        global <- result$common_component$global
        expect_identical(dim(global$G_hat), c(12L, kr, kc))
        expect_identical(dim(global$G_tilde), c(12L, kr, kc))
        for (m in seq_along(fixture$Y)) {
          expect_identical(dim(global$Q_tilde[[m]]), c(dim(fixture$Y[[m]])[2L], kr))
          expect_identical(dim(global$J_tilde[[m]]), c(dim(fixture$Y[[m]])[3L], kc))
          expect_equal(.common_test_product(global$Q_tilde[[m]], global$G_tilde, global$J_tilde[[m]]),
                       fixture$global[[m]], tolerance = 1e-10)
        }
      }
    }
  }
})
