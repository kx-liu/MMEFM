# Exact common-only fixtures with unequal group dimensions and centred loading columns.
.common_projection_fixture <- function(kr = 1L, kc = 1L) {
  TT <- 8L
  p <- c(4L, 5L, 6L)
  q <- c(5L, 6L, 4L)
  temporal <- unname(contr.helmert(TT))
  temporal <- sweep(temporal, 2L, sqrt(colSums(temporal^2)), "/")
  check_Y <- Q <- J <- vector("list", 3L)
  names(check_Y) <- names(Q) <- names(J) <- c("a", "b", "c")
  for (m in seq_len(3L)) {
    Q[[m]] <- unname(contr.helmert(p[m]))[, seq_len(kr), drop = FALSE]
    J[[m]] <- unname(contr.helmert(q[m]))[, seq_len(kc), drop = FALSE]
    Q[[m]] <- (m + 1) * sweep(Q[[m]], 2L, sqrt(colSums(Q[[m]]^2)), "/")
    J[[m]] <- (m + 2) * sweep(J[[m]], 2L, sqrt(colSums(J[[m]]^2)), "/")
    check_Y[[m]] <- array(0, c(TT, p[m], q[m]))
    for (t in seq_len(TT)) {
      G <- matrix(c(5, 3, 2, 1)[seq_len(kr * kc)] * temporal[t, seq_len(kr * kc)], kr, kc)
      check_Y[[m]][t, , ] <- Q[[m]] %*% G %*% t(J[[m]])
    }
  }
  list(check_Y = check_Y, Q = Q, J = J)
}

test_that("global common PSDs match the explicit time-sum formula on both sides", {
  for (TT in c(1L, 3L)) {
    check_Y <- list(a = array(sin(seq_len(TT * 3L * 4L)), c(TT, 3L, 4L)), b = array(cos(seq_len(TT * 4L * 3L)), c(TT, 4L, 3L)), c = array(sin(2 * seq_len(TT * 2L * 5L)), c(TT, 2L, 5L)))
    for (side in c("row", "column")) {
      directions <- vector("list", length(check_Y))
      for (m in seq_along(check_Y)) {
        d <- dim(check_Y[[m]])[if (side == "row") 3L else 2L]
        directions[[m]] <- matrix(seq_len(d) / sqrt(sum(seq_len(d)^2)), d, 1L)
      }
      psd <- .common_global_psd(check_Y, directions, side)
      expect_identical(names(psd), names(check_Y))
      for (m in seq_along(check_Y)) {
        d_m <- dim(check_Y[[m]])[if (side == "row") 2L else 3L]
        expected <- matrix(0, d_m, d_m)
        for (n in setdiff(seq_along(check_Y), m)) {
          d_n <- dim(check_Y[[n]])[if (side == "row") 2L else 3L]
          covariance_sum <- matrix(0, d_m, d_n)
          for (t in seq_len(TT)) {
            Ym <- matrix(check_Y[[m]][t, , , drop = FALSE], dim(check_Y[[m]])[2L], dim(check_Y[[m]])[3L])
            Yn <- matrix(check_Y[[n]][t, , , drop = FALSE], dim(check_Y[[n]])[2L], dim(check_Y[[n]])[3L])
            if (side == "row") {
              covariance_sum <- covariance_sum + Ym %*% directions[[m]] %*% t(directions[[n]]) %*% t(Yn)
            } else {
              covariance_sum <- covariance_sum + t(Ym) %*% directions[[m]] %*% t(directions[[n]]) %*% Yn
            }
          }
          expected <- expected + tcrossprod(covariance_sum) / (TT^2 * (length(check_Y) - 1L))
        }
        expect_identical(dim(psd[[m]]), c(d_m, d_m))
        expect_identical(psd[[m]], t(psd[[m]]))
        expect_equal(psd[[m]], expected, tolerance = 1e-12)
      }
    }
  }
})

test_that("direction search uses updates from preceding groups", {
  # Group 2 prefers candidate 1 against the initial group-1 direction, but candidate 2 after group 1 updates.
  check_Y <- list(a = array(0, c(2L, 2L, 2L)), b = array(0, c(2L, 2L, 2L)))
  check_Y$a[, 1L, 1L] <- c(1, 0)
  check_Y$a[, 1L, 2L] <- c(0, 2)
  check_Y$b[, 1L, 1L] <- c(1, 1)
  check_Y$b[, 1L, 2L] <- c(0, 3)
  pools <- list(list(matrix(c(1, 0), 2L, 1L), matrix(c(0, 1), 2L, 1L)), list(matrix(c(1, 0), 2L, 1L), matrix(c(0, 1), 2L, 1L)))
  for (side in c("row", "column")) {
    Y <- if (side == "row") check_Y else lapply(check_Y, aperm, perm = c(1L, 3L, 2L))
    directions <- .select_common_directions(Y, pools, side)
    expect_identical(names(directions), names(check_Y))
    for (m in seq_along(Y)) expect_identical(directions[[m]], pools[[m]][[2L]])
  }
  # An exact tie uses the ordinary first-maximum rule.
  zero <- lapply(check_Y, function(Y) array(0, dim(Y)))
  directions <- .select_common_directions(zero, pools, "row")
  for (m in seq_along(zero)) expect_identical(directions[[m]], pools[[m]][[1L]])
})

test_that("rank-one Algorithm 1 recovers spaces reproducibly with caller-controlled RNG", {
  fixture <- .common_projection_fixture()
  set.seed(2026L)
  first <- .initial_common_loadings(fixture$check_Y, kr = 1L, kc = 1L, K0 = 8L)
  set.seed(2026L)
  second <- .initial_common_loadings(fixture$check_Y, kr = 1L, kc = 1L, K0 = 8L)
  expect_identical(names(first), c("Q_hat", "J_hat", "iterations", "converged"))
  expect_identical(names(first$Q_hat), names(fixture$check_Y))
  expect_identical(names(first$J_hat), names(fixture$check_Y))
  expect_type(first$iterations, "integer")
  expect_true(first$iterations >= 1L && first$iterations <= 20L)
  expect_true(first$converged)
  expect_identical(first$iterations, second$iterations)
  expect_identical(first$converged, second$converged)
  for (m in seq_along(fixture$check_Y)) {
    expect_identical(dim(first$Q_hat[[m]]), c(dim(fixture$check_Y[[m]])[2L], 1L))
    expect_identical(dim(first$J_hat[[m]]), c(dim(fixture$check_Y[[m]])[3L], 1L))
    expect_equal(tcrossprod(first$Q_hat[[m]]), tcrossprod(fixture$Q[[m]]) / sum(fixture$Q[[m]]^2), tolerance = 1e-12)
    expect_equal(tcrossprod(first$J_hat[[m]]), tcrossprod(fixture$J[[m]]) / sum(fixture$J[[m]]^2), tolerance = 1e-12)
    expect_identical(tcrossprod(first$Q_hat[[m]]), tcrossprod(second$Q_hat[[m]]))
    expect_identical(tcrossprod(first$J_hat[[m]]), tcrossprod(second$J_hat[[m]]))
  }

  # Only the M*K0 half-panel draws advance the RNG; the initializer never resets it.
  set.seed(907L)
  for (m in seq_along(fixture$check_Y)) {
    q <- dim(fixture$check_Y[[m]])[3L]
    for (k in seq_len(8L)) sample(q, floor(q / 2))
  }
  expected_seed <- .Random.seed
  set.seed(907L)
  invisible(.initial_common_loadings(fixture$check_Y, 1L, 1L, K0 = 8L))
  expect_identical(.Random.seed, expected_seed)
  expect_false("seed" %in% names(formals(.initial_common_loadings)))
})

test_that("a nondegenerate two-by-two common fixture recovers both loading spaces", {
  fixture <- .common_projection_fixture(kr = 2L, kc = 2L)
  set.seed(2026L)
  fit <- .initial_common_loadings(fixture$check_Y, kr = 2L, kc = 2L, K0 = 8L)
  expect_true(fit$converged)
  for (m in seq_along(fixture$check_Y)) {
    Q <- fixture$Q[[m]] / (m + 1)
    J <- fixture$J[[m]] / (m + 2)
    expect_identical(dim(fit$Q_hat[[m]]), c(nrow(Q), 2L))
    expect_identical(dim(fit$J_hat[[m]]), c(nrow(J), 2L))
    expect_equal(tcrossprod(fit$Q_hat[[m]]), tcrossprod(Q), tolerance = 1e-12)
    expect_equal(tcrossprod(fit$J_hat[[m]]), tcrossprod(J), tolerance = 1e-12)
    expect_equal(crossprod(fit$Q_hat[[m]]), diag(2L), tolerance = 1e-12)
    expect_equal(crossprod(fit$J_hat[[m]]), diag(2L), tolerance = 1e-12)
  }
})

test_that("principal-angle distance is invariant to signs and basis rotations", {
  A <- diag(3L)[, 1:2, drop = FALSE]
  rotation <- matrix(c(0, 1, -1, 0), 2L, 2L)
  expect_equal(.subspace_distance(A, -A), 0, tolerance = 1e-12)
  expect_equal(.subspace_distance(A, A %*% rotation), 0, tolerance = 1e-12)
  expect_equal(.subspace_distance(A[, 1L, drop = FALSE], diag(3L)[, 3L, drop = FALSE]), 1, tolerance = 1e-12)
})
