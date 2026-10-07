# Shared exact additive fixture; global/local loading overlap exercises the projected Gram systems.
.main_effect_fixture <- function() {
  TT <- 12L
  p <- c(4L, 5L, 6L)
  q <- c(5L, 4L, 6L)
  temporal <- unname(contr.helmert(TT))
  temporal <- sweep(temporal, 2L, sqrt(colSums(temporal^2) / TT), "/")
  alpha <- temporal[, 1L, drop = FALSE]
  beta <- temporal[, 2L, drop = FALSE]
  mu <- matrix(2, TT, 1L)
  Xt <- A1 <- A2 <- B1 <- B2 <- local_mu <- vector("list", 3L)
  names(Xt) <- names(A1) <- names(A2) <- names(B1) <- names(B2) <- names(local_mu) <- c("a", "b", "c")
  for (m in seq_len(3L)) {
    H_A <- unname(contr.helmert(p[m]))
    H_B <- unname(contr.helmert(q[m]))
    A1[[m]] <- (m + 1) * H_A[, 1L, drop = FALSE]
    B1[[m]] <- (m + 2) * H_B[, 1L, drop = FALSE]
    A2[[m]] <- H_A[, 2L, drop = FALSE] + 0.35 * H_A[, 1L, drop = FALSE]
    B2[[m]] <- H_B[, 2L, drop = FALSE] - 0.25 * H_B[, 1L, drop = FALSE]
    local_mu[[m]] <- (m - 2) * temporal[, 9L, drop = FALSE]
    a <- alpha %*% t(A1[[m]]) + temporal[, m + 2L, drop = FALSE] %*% t(A2[[m]])
    b <- beta %*% t(B1[[m]]) + temporal[, m + 5L, drop = FALSE] %*% t(B2[[m]])
    Xt[[m]] <- array(0, c(TT, p[m], q[m]))
    for (t in seq_len(TT)) {
      Xt[[m]][t, , ] <- matrix(mu[t, 1L] + local_mu[[m]][t, 1L], p[m], q[m]) + outer(a[t, ], rep(1, q[m])) + outer(rep(1, p[m]), b[t, ])
    }
  }
  rank <- .validate_rank(list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3L), l2 = rep(1L, 3L), kr = 1L, kc = 1L, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L)), .validate_Xt(Xt))
  list(Xt = Xt, rank = rank, A1 = A1, A2 = A2, B1 = B1, B2 = B2, mu = mu, local_mu = local_mu)
}

test_that("direct moments satisfy centering, residual, and singleton-time identities", {
  for (TT in c(1L, 2L)) {
    Xt <- list(a = array(sin(seq_len(TT * 3L * 4L)), c(TT, 3L, 4L)), b = array(cos(seq_len(TT * 4L * 3L)), c(TT, 4L, 3L)))
    moments <- .main_effect_moments(Xt)
    expect_identical(names(moments), c("direct", "check_Y"))
    for (m in seq_along(Xt)) {
      p <- dim(Xt[[m]])[2L]
      q <- dim(Xt[[m]])[3L]
      direct <- moments$direct
      expect_identical(dim(direct$c_hat[[m]]), c(TT, 1L))
      expect_identical(dim(direct$a_hat[[m]]), c(TT, p))
      expect_identical(dim(direct$b_hat[[m]]), c(TT, q))
      expect_identical(dim(moments$check_Y[[m]]), dim(Xt[[m]]))
      expect_equal(rowMeans(direct$a_hat[[m]]), rep(0, TT), tolerance = 1e-12)
      expect_equal(rowMeans(direct$b_hat[[m]]), rep(0, TT), tolerance = 1e-12)
      for (t in seq_len(TT)) {
        Y <- matrix(Xt[[m]][t, , , drop = FALSE], p, q)
        residual <- matrix(moments$check_Y[[m]][t, , , drop = FALSE], p, q)
        expect_equal(direct$c_hat[[m]][t, 1L], mean(Y), tolerance = 1e-12)
        expect_equal(direct$a_hat[[m]][t, ], rowMeans(Y) - mean(Y), tolerance = 1e-12)
        expect_equal(direct$b_hat[[m]][t, ], colMeans(Y) - mean(Y), tolerance = 1e-12)
        expected <- Y - matrix(direct$c_hat[[m]][t, 1L], p, q) - outer(direct$a_hat[[m]][t, ], rep(1, q)) - outer(rep(1, p), direct$b_hat[[m]][t, ])
        expect_equal(residual, expected, tolerance = 1e-12)
        expect_equal(rowMeans(residual), rep(0, p), tolerance = 1e-12)
        expect_equal(colMeans(residual), rep(0, q), tolerance = 1e-12)
      }
    }
  }
})

test_that("main-effect PSD matrices use the paper covariance normalizations and projections", {
  fixture <- .main_effect_fixture()
  direct <- .main_effect_moments(fixture$Xt)$direct
  TT <- dim(fixture$Xt[[1]])[1L]
  psd <- .main_effect_global_psd(direct$a_hat, direct$b_hat, TT)
  local <- .main_effect_local_psd(direct$a_hat, direct$b_hat, lapply(fixture$A1, function(A) A / sqrt(sum(A^2))), lapply(fixture$B1, function(B) B / sqrt(sum(B^2))), TT)
  for (m in seq_along(fixture$Xt)) {
    expected_A <- matrix(0, nrow(fixture$A1[[m]]), nrow(fixture$A1[[m]]))
    expected_B <- matrix(0, nrow(fixture$B1[[m]]), nrow(fixture$B1[[m]]))
    for (n in setdiff(seq_along(fixture$Xt), m)) {
      expected_A <- expected_A + sum(fixture$A1[[n]]^2) * tcrossprod(fixture$A1[[m]]) / 2
      expected_B <- expected_B + sum(fixture$B1[[n]]^2) * tcrossprod(fixture$B1[[m]]) / 2
    }
    expect_equal(psd$row[[m]], expected_A, tolerance = 1e-12)
    expect_equal(psd$column[[m]], expected_B, tolerance = 1e-12)
    P_A <- diag(nrow(fixture$A1[[m]])) - tcrossprod(fixture$A1[[m]]) / sum(fixture$A1[[m]]^2)
    P_B <- diag(nrow(fixture$B1[[m]])) - tcrossprod(fixture$B1[[m]]) / sum(fixture$B1[[m]]^2)
    expect_equal(local$row[[m]], sum(crossprod(fixture$A2[[m]], P_A %*% fixture$A2[[m]])) * tcrossprod(fixture$A2[[m]]), tolerance = 1e-12)
    expect_equal(local$column[[m]], sum(crossprod(fixture$B2[[m]], P_B %*% fixture$B2[[m]])) * tcrossprod(fixture$B2[[m]]), tolerance = 1e-12)
  }
  expect_identical(names(psd$row), names(fixture$Xt))
  expect_identical(names(psd$column), names(fixture$Xt))
  expect_identical(names(local$row), names(fixture$Xt))
  expect_identical(names(local$column), names(fixture$Xt))
})

test_that("the exact additive model recovers loading spaces, normalized scores, and the paper refit", {
  fixture <- .main_effect_fixture()
  fit <- .estimate_main_effect(fixture$Xt, fixture$rank)
  global <- fit$main_effect$global
  local <- fit$main_effect$local
  direct <- fit$main_effect$direct
  TT <- dim(fixture$Xt[[1]])[1L]
  expect_identical(names(fit), c("main_effect", "check_Y"))
  expect_identical(names(fit$main_effect), c("direct", "global", "local"))
  expect_identical(names(direct), c("c_hat", "a_hat", "b_hat"))
  expect_identical(names(global), c("mu_hat", "A1_hat", "alpha_hat", "A1_tilde", "B1_hat", "beta_hat", "B1_tilde"))
  expect_identical(names(local), c("mu_hat", "A2_hat", "alpha_hat", "B2_hat", "beta_hat"))
  expect_identical(dim(global$mu_hat), c(TT, 1L))
  expect_identical(dim(global$alpha_hat), c(TT, 1L))
  expect_identical(dim(global$beta_hat), c(TT, 1L))
  expect_equal(crossprod(global$alpha_hat) / TT, diag(1L), tolerance = 1e-12)
  expect_equal(crossprod(global$beta_hat) / TT, diag(1L), tolerance = 1e-12)
  expect_equal(global$mu_hat, fixture$mu, tolerance = 1e-12)
  for (m in seq_along(fixture$Xt)) {
    p <- dim(fixture$Xt[[m]])[2L]
    q <- dim(fixture$Xt[[m]])[3L]
    for (A in c("A1_hat", "A1_tilde")) expect_equal(tcrossprod(global[[A]][[m]]) / sum(global[[A]][[m]]^2), tcrossprod(fixture$A1[[m]]) / sum(fixture$A1[[m]]^2), tolerance = 1e-12)
    for (B in c("B1_hat", "B1_tilde")) expect_equal(tcrossprod(global[[B]][[m]]) / sum(global[[B]][[m]]^2), tcrossprod(fixture$B1[[m]]) / sum(fixture$B1[[m]]^2), tolerance = 1e-12)
    expect_equal(tcrossprod(local$A2_hat[[m]]), tcrossprod(fixture$A2[[m]]) / sum(fixture$A2[[m]]^2), tolerance = 1e-12)
    expect_equal(tcrossprod(local$B2_hat[[m]]), tcrossprod(fixture$B2[[m]]) / sum(fixture$B2[[m]]^2), tolerance = 1e-12)
    for (A in c("A1_hat", "A1_tilde")) expect_identical(dim(global[[A]][[m]]), c(p, 1L))
    for (B in c("B1_hat", "B1_tilde")) expect_identical(dim(global[[B]][[m]]), c(q, 1L))
    expect_identical(dim(local$A2_hat[[m]]), c(p, 1L))
    expect_identical(dim(local$B2_hat[[m]]), c(q, 1L))
    expect_identical(dim(local$alpha_hat[[m]]), c(TT, 1L))
    expect_identical(dim(local$beta_hat[[m]]), c(TT, 1L))
    expect_identical(dim(local$mu_hat[[m]]), c(TT, 1L))
    expect_equal(local$mu_hat[[m]], fixture$local_mu[[m]], tolerance = 1e-12)
    local_a <- local$alpha_hat[[m]] %*% t(local$A2_hat[[m]])
    local_b <- local$beta_hat[[m]] %*% t(local$B2_hat[[m]])
    y_hat <- direct$a_hat[[m]] - local_a
    w_hat <- direct$b_hat[[m]] - local_b
    expect_equal(global$A1_tilde[[m]], crossprod(y_hat, global$alpha_hat) / TT, tolerance = 1e-12)
    expect_equal(global$B1_tilde[[m]], crossprod(w_hat, global$beta_hat) / TT, tolerance = 1e-12)
    a <- global$alpha_hat %*% t(global$A1_tilde[[m]]) + local_a
    b <- global$beta_hat %*% t(global$B1_tilde[[m]]) + local_b
    for (t in seq_len(TT)) {
      reconstructed <- matrix(global$mu_hat[t, 1L] + local$mu_hat[[m]][t, 1L], p, q) + outer(a[t, ], rep(1, q)) + outer(rep(1, p), b[t, ])
      expect_equal(reconstructed, matrix(fixture$Xt[[m]][t, , , drop = FALSE], p, q), tolerance = 1e-12)
    }
    expect_equal(fit$check_Y[[m]], array(0, dim(fixture$Xt[[m]])), tolerance = 1e-12)
  }
  for (group_list in c(direct, global[c("A1_hat", "A1_tilde", "B1_hat", "B1_tilde")], local, list(fit$check_Y))) expect_identical(names(group_list), names(fixture$Xt))
  # Common ranks are irrelevant to this internal main-effect calculation.
  expect_identical(.estimate_main_effect(fixture$Xt, fixture$rank[c("r1", "l1", "r2", "l2")]), fit)
})

test_that("the optional ridge uses explicit regularized Gram systems", {
  fixture <- .main_effect_fixture()
  lambda <- 0.1
  fit <- .estimate_main_effect(fixture$Xt, fixture$rank, lambda)
  global <- fit$main_effect$global
  local <- fit$main_effect$local
  direct <- fit$main_effect$direct
  expect_true(all(is.finite(unlist(fit))))
  expect_identical(dim(global$alpha_hat), c(12L, 1L))
  expect_identical(dim(global$beta_hat), c(12L, 1L))
  for (m in seq_along(fixture$Xt)) {
    p <- dim(fixture$Xt[[m]])[2L]
    q <- dim(fixture$Xt[[m]])[3L]
    P_A <- diag(p) - tcrossprod(global$A1_hat[[m]])
    P_B <- diag(q) - tcrossprod(global$B1_hat[[m]])
    A2 <- local$A2_hat[[m]]
    B2 <- local$B2_hat[[m]]
    expect_equal(local$alpha_hat[[m]], t(solve(crossprod(A2, P_A %*% A2) + lambda * diag(1L), t(A2) %*% P_A %*% t(direct$a_hat[[m]]))), tolerance = 1e-12)
    expect_equal(local$beta_hat[[m]], t(solve(crossprod(B2, P_B %*% B2) + lambda * diag(1L), t(B2) %*% P_B %*% t(direct$b_hat[[m]]))), tolerance = 1e-12)
    y_hat <- direct$a_hat[[m]] - local$alpha_hat[[m]] %*% t(A2)
    w_hat <- direct$b_hat[[m]] - local$beta_hat[[m]] %*% t(B2)
    expect_equal(global$A1_tilde[[m]], crossprod(y_hat, global$alpha_hat) %*% solve(crossprod(global$alpha_hat) + lambda * diag(1L)), tolerance = 1e-12)
    expect_equal(global$B1_tilde[[m]], crossprod(w_hat, global$beta_hat) %*% solve(crossprod(global$beta_hat) + lambda * diag(1L)), tolerance = 1e-12)
    expect_identical(dim(global$A1_tilde[[m]]), c(p, 1L))
    expect_identical(dim(global$B1_tilde[[m]]), c(q, 1L))
    expect_identical(dim(local$alpha_hat[[m]]), c(12L, 1L))
    expect_identical(dim(local$beta_hat[[m]]), c(12L, 1L))
  }
  for (lambda in c(-1, NA_real_, Inf)) expect_error(.estimate_main_effect(fixture$Xt, fixture$rank, lambda), "finite nonnegative")
})

test_that("numerical primitives symmetrize eigenproblems and expose singular solves", {
  S <- matrix(c(3, 1, -1, 2), 2L)
  V <- .leading_eigenvectors(S, 1L)
  expect_identical(dim(V), c(2L, 1L))
  expect_equal(tcrossprod(V), diag(c(1, 0)), tolerance = 1e-12)
  A <- matrix(1, 2L, 2L)
  B <- diag(2L)
  expect_error(.solve_gram(A, B), "singular")
  expect_equal(.solve_gram(A, B, lambda = 0.1), solve(A + 0.1 * diag(2L), B), tolerance = 1e-12)
})
