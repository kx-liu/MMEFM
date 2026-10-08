.rank_selection_fixture <- function() {
  set.seed(917L)
  list(a = array(rnorm(12L * 2L * 2L), c(12L, 2L, 2L)),
       b = array(rnorm(12L * 2L * 2L), c(12L, 2L, 2L)),
       c = array(rnorm(12L * 2L * 2L), c(12L, 2L, 2L)))
}

# Minimum dimensions leave the true elbow as the last available ratio; no numerical perturbation is needed.
.rank_common_fixture <- function(kr = 1L, kc = 1L) {
  TT <- 6L
  p <- kr + c(1L, 2L, 3L)
  q <- kc + c(1L, 2L, 3L)
  temporal <- unname(contr.helmert(TT))
  temporal <- sweep(temporal, 2L, sqrt(colSums(temporal^2)), "/")
  Q <- J <- Y <- vector("list", 3L)
  names(Q) <- names(J) <- names(Y) <- c("a", "b", "c")
  for (m in seq_len(3L)) {
    Q[[m]] <- unname(contr.helmert(p[m]))[, seq_len(kr), drop = FALSE]
    J[[m]] <- unname(contr.helmert(q[m]))[, seq_len(kc), drop = FALSE]
    Q[[m]] <- sweep(Q[[m]], 2L, sqrt(colSums(Q[[m]]^2)), "/")
    J[[m]] <- sweep(J[[m]], 2L, sqrt(colSums(J[[m]]^2)), "/")
    Y[[m]] <- array(0, c(TT, p[m], q[m]))
    for (t in seq_len(TT)) {
      Gt <- matrix(c(5, 2)[seq_len(kr * kc)] * temporal[t, seq_len(kr * kc)], kr, kc)
      Y[[m]][t, , ] <- Q[[m]] %*% Gt %*% t(J[[m]])
    }
  }
  list(Y = Y, Q = Q, J = J)
}

test_that("local paper ratios use every consecutive eigenvalue without perturbation", {
  spectrum <- c(100, 80, 60, 40, 30, 0.01, 0.009)
  fit <- .paper_eigen_ratio(diag(spectrum))
  expect_identical(fit$rank, 5L)
  expect_equal(fit$ratios, spectrum[-1L] / spectrum[-length(spectrum)], tolerance = 1e-15)
  expect_length(fit$ratios, length(spectrum) - 1L)
  scaled <- .paper_eigen_ratio(8 * diag(spectrum))
  expect_identical(scaled$rank, fit$rank)
  expect_equal(scaled$ratios, fit$ratios, tolerance = 1e-15)

  zeros <- .paper_eigen_ratio(diag(c(4, 1, 0, 0)))
  expect_identical(zeros$rank, 2L)
  expect_identical(zeros$ratios, c(0.25, 0, Inf))
  all_zero <- .paper_eigen_ratio(matrix(0, 3L, 3L))
  expect_identical(all_zero$rank, 1L)
  expect_identical(all_zero$ratios, c(Inf, Inf))
  expect_identical(all_zero$pooled_eigenvalues, rep(0, 3L))
  roundoff <- .paper_eigen_ratio(diag(c(4, 1, -.Machine$double.eps)))
  expect_equal(roundoff$eigenvalues[[1L]], c(4, 1, 0))
  expect_error(.paper_eigen_ratio(diag(c(4, 1, -1e-7))), "materially")
  asymmetric <- diag(c(4, 2, 1))
  asymmetric[1L, 2L] <- 1
  asymmetric[2L, 1L] <- -1
  expect_equal(.paper_eigen_ratio(asymmetric)$pooled_eigenvalues, c(4, 2, 1))
  expect_error(.paper_eigen_ratio(list()), "nonempty")
  expect_error(.paper_eigen_ratio(matrix(1, 1L, 1L)), "at least two")
})

test_that("global ranks use indexwise spectral maxima over the shared range", {
  matrices <- list(a = diag(c(101, 100, 1, 0.9)), b = diag(c(100, 50, 49, 0.1, 0.01)))
  fit <- .paper_eigen_ratio(matrices)
  expect_identical(names(fit$eigenvalues), names(matrices))
  expect_identical(names(fit$eigenvectors), names(matrices))
  expect_equal(fit$pooled_eigenvalues, c(101, 100, 49, 0.9))
  for (j in seq_len(4L)) expect_equal(fit$pooled_eigenvalues[j], max(fit$eigenvalues$a[j], fit$eigenvalues$b[j]))
  expect_equal(fit$ratios, c(100 / 101, 49 / 100, 0.9 / 49))
  expect_identical(fit$rank, 3L)
  expect_identical(.paper_eigen_ratio(matrices$a)$rank, 2L)
  expect_length(fit$ratios, 3L)
})

test_that("automatic Algorithm 1 recovers exact ranks and preserves reproducibility", {
  for (ranks in list(c(1L, 1L), c(2L, 1L), c(1L, 2L))) {
    fixture <- .rank_common_fixture(ranks[1L], ranks[2L])
    set.seed(2026L)
    first <- .initial_common_loadings(fixture$Y, K0 = 8L)
    set.seed(2026L)
    second <- .initial_common_loadings(fixture$Y, K0 = 8L)
    expect_identical(first, second)
    expect_identical(c(first$kr, first$kc), ranks)
    expect_true(first$converged)
    expect_identical(first$row_rank_fit$rank, first$kr)
    expect_identical(first$column_rank_fit$rank, first$kc)
    expect_identical(names(first), c("Q_hat", "J_hat", "iterations", "converged", "kr", "kc", "row_rank_fit", "column_rank_fit"))
    for (m in seq_along(fixture$Y)) {
      expect_identical(dim(first$Q_hat[[m]]), c(nrow(fixture$Q[[m]]), ranks[1L]))
      expect_identical(dim(first$J_hat[[m]]), c(nrow(fixture$J[[m]]), ranks[2L]))
      expect_equal(tcrossprod(first$Q_hat[[m]]), tcrossprod(fixture$Q[[m]]), tolerance = 1e-12)
      expect_equal(tcrossprod(first$J_hat[[m]]), tcrossprod(fixture$J[[m]]), tolerance = 1e-12)
    }
  }
  expect_error(.initial_common_loadings(fixture$Y, kr = 1L), "both")
  expect_error(.initial_common_loadings(fixture$Y, kc = 1L), "both")
  expect_false("seed" %in% names(formals(.initial_common_loadings)))
  set.seed(2026L)
  before <- .Random.seed
  invisible(.initial_common_loadings(fixture$Y, K0 = 8L))
  expect_false(identical(.Random.seed, before))
})

test_that("a rank change prevents projection convergence even with a large space tolerance", {
  set.seed(1L)
  Y <- list(a = array(rnorm(8L * 4L * 5L), c(8L, 4L, 5L)),
            b = array(rnorm(8L * 5L * 4L), c(8L, 5L, 4L)),
            c = array(rnorm(8L * 6L * 6L), c(8L, 6L, 6L)))
  set.seed(2026L)
  fit <- .initial_common_loadings(Y, K0 = 4L, max_iter = 1L, tol = 1e9)
  expect_false(fit$converged)
  expect_identical(c(fit$kr, fit$kc), c(1L, 1L))
  expect_identical(fit$iterations, 1L)
})

test_that("public rank objects and diagnostics match the canonical contracts", {
  Xt <- .rank_selection_fixture()
  result <- select_MMEFM_rank(Xt, K0 = 8L)
  components <- c("r1", "l1", "r2", "l2", "kr", "kc", "kr_m", "kc_m")
  expect_s3_class(result, "mmefm_rank")
  expect_identical(names(result), c(components, "diagnostics"))
  expect_identical(.validate_rank(result, .validate_Xt(Xt)), unclass(result)[components])
  for (field in c("r1", "l1", "kr", "kc")) {
    expect_type(result[[field]], "integer")
    expect_length(result[[field]], 1L)
  }
  for (field in c("r2", "l2", "kr_m", "kc_m")) {
    expect_type(result[[field]], "integer")
    expect_length(result[[field]], length(Xt))
    expect_identical(names(result[[field]]), names(Xt))
  }
  diagnostics <- result$diagnostics
  expect_identical(names(diagnostics), components)
  fit_names <- c("eigenvalues", "pooled_eigenvalues", "ratios")
  for (field in c("r1", "l1", "kr", "kc")) {
    expect_identical(names(diagnostics[[field]]), fit_names)
    expect_identical(names(diagnostics[[field]]$eigenvalues), names(Xt))
    expect_identical(as.integer(which.min(diagnostics[[field]]$ratios)), result[[field]])
  }
  for (field in c("r2", "l2", "kr_m", "kc_m")) {
    expect_identical(names(diagnostics[[field]]), names(Xt))
    for (m in seq_along(Xt)) {
      expect_identical(names(diagnostics[[field]][[m]]), c("eigenvalues", "ratios"))
      expect_identical(as.integer(which.min(diagnostics[[field]][[m]]$ratios)), unname(result[[field]][m]))
      expect_length(diagnostics[[field]][[m]]$ratios, length(diagnostics[[field]][[m]]$eigenvalues) - 1L)
    }
  }
  expect_null(result$source)
  unnamed <- select_MMEFM_rank(unname(Xt), K0 = 8L)
  for (field in c("r2", "l2", "kr_m", "kc_m")) expect_null(names(unnamed[[field]]))
})

test_that("the public selector restores caller RNG and ignores its initial state", {
  Xt <- .rank_selection_fixture()
  set.seed(71L)
  before <- .Random.seed
  kind <- RNGkind()
  first <- select_MMEFM_rank(Xt, K0 = 8L, seed = 2026L)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  set.seed(903L)
  before <- .Random.seed
  second <- select_MMEFM_rank(Xt, K0 = 8L, seed = 2026L)
  expect_identical(.Random.seed, before)
  expect_identical(first, second)
  rm(".Random.seed", envir = .GlobalEnv)
  third <- select_MMEFM_rank(Xt, K0 = 8L, seed = 2026L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(first, third)
})

test_that("the public selector restores RNG after an error in its stochastic section", {
  # Direct moments cancel exactly; only the stochastic common cross-group Gram overflows.
  Y <- array(1e100 * rep(c(1, -1, -1, 1), each = 2L), c(2L, 2L, 2L))
  Xt <- list(a = Y, b = Y)
  set.seed(71L)
  before <- .Random.seed
  expect_error(select_MMEFM_rank(Xt), "infinite or missing")
  expect_identical(.Random.seed, before)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_error(select_MMEFM_rank(Xt), "infinite or missing")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("public rank controls are validated without rounding", {
  Xt <- .rank_selection_fixture()
  for (value in list(0, 1.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    expect_error(select_MMEFM_rank(Xt, K0 = value), "K0")
    expect_error(select_MMEFM_rank(Xt, max_iter = value), "max_iter")
  }
  for (value in list(0, -1, Inf, NA_real_)) expect_error(select_MMEFM_rank(Xt, tol = value), "tol")
  for (value in list(-1, 1.5, Inf, NA_real_, .Machine$integer.max + 1)) expect_error(select_MMEFM_rank(Xt, seed = value), "seed")
  for (value in list(1, NA, c(TRUE, FALSE))) expect_error(select_MMEFM_rank(Xt, verbose = value), "verbose")
  expect_error(select_MMEFM_rank(Xt[1L]), "at least two")
  expect_identical(select_MMEFM_rank(Xt, K0 = 2, max_iter = 1, seed = 0),
                   select_MMEFM_rank(Xt, K0 = 2L, max_iter = 1L, seed = 0L))
})
