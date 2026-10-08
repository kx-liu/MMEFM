.fit_fixture <- function(p = c(4L, 5L, 6L), q = c(5L, 6L, 4L), global = c(1, -1, 0), local = c(1, 1, -2)) {
  H <- matrix(1, 1L, 1L)
  for (i in seq_len(4L)) H <- rbind(cbind(H, H), cbind(H, -H))
  Xt <- global_main <- local_main <- global_common <- local_common <- vector("list", 3L)
  names(Xt) <- names(global_main) <- names(local_main) <- names(global_common) <- names(local_common) <- c("a", "b", "c")
  for (m in seq_len(3L)) {
    a <- matrix(c(global, rep(0, p[m] - length(global))), p[m], 1L)
    b <- matrix(c(global, rep(0, q[m] - length(global))), q[m], 1L)
    r <- matrix(c(local, rep(0, p[m] - length(local))), p[m], 1L)
    c <- matrix(c(local, rep(0, q[m] - length(local))), q[m], 1L)
    global_main[[m]] <- local_main[[m]] <- global_common[[m]] <- local_common[[m]] <- array(0, c(16L, p[m], q[m]))
    for (t in seq_len(16L)) {
      global_main[[m]][t, , ] <- matrix(H[t, 13L], p[m], q[m]) +
        (m + 1) * 3 * H[t, 2L] * a %*% matrix(1, 1L, q[m]) +
        (m + 2) * 2 * H[t, 6L] * matrix(1, p[m], 1L) %*% t(b)
      local_main[[m]][t, , ] <- matrix((m - 2) * H[t, 14L], p[m], q[m]) +
        H[t, 2L + m] * r %*% matrix(1, 1L, q[m]) + H[t, 6L + m] * matrix(1, p[m], 1L) %*% t(c)
      global_common[[m]][t, , ] <- (m + 1) * (m + 2) * 5 * H[t, 10L] * a %*% t(b)
      local_common[[m]][t, , ] <- H[t, 10L + m] * r %*% t(c)
    }
    Xt[[m]] <- global_main[[m]] + local_main[[m]] + global_common[[m]] + local_common[[m]]
  }
  list(Xt = Xt, global_main = global_main, local_main = local_main, global_common = global_common, local_common = local_common,
       rank = list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3L), l2 = rep(1L, 3L), kr = 1L, kc = 1L, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L)))
}

test_that("complete fits reconstruct all identified components", {
  x <- .fit_fixture()
  for (alignment in c("VanLoan", "Procrustes")) {
    for (refit in c("VanLoan", "ALS")) {
      expect_warning(fit <- est_MMEFM(x$Xt, x$rank, K0 = 8L, alignment_method = alignment, refit_method = refit), NA)
      expect_s3_class(fit, "mmefm_fit")
      expect_identical(names(fit), c("call", "dimensions", "rank", "main_effect", "common_component", "check_Y", "convergence"))
      expect_s3_class(fit$rank, "mmefm_rank")
      expect_null(fit$rank$diagnostics)
      expect_identical(names(fit$dimensions$p), names(x$Xt))
      expect_identical(names(fit$dimensions$q), names(x$Xt))
      for (name in c("r2", "l2", "kr_m", "kc_m")) expect_identical(names(fit$rank[[name]]), names(x$Xt))
      for (component in c("global_main", "local_main", "global_common", "local_common", "main", "all")) {
        values <- fitted(fit, component)
        truth <- switch(component, main = Map(`+`, x$global_main, x$local_main), all = x$Xt, x[[component]])
        expect_identical(names(values), names(x$Xt))
        for (m in seq_along(values)) {
          expect_identical(dim(values[[m]]), dim(x$Xt[[m]]))
          expect_equal(values[[m]], truth[[m]], tolerance = 1e-10)
        }
      }
      errors <- residuals(fit)
      expect_identical(names(errors), names(x$Xt))
      for (m in seq_along(x$Xt)) {
        observed <- array(0, dim(x$Xt[[m]]))
        for (t in seq_len(dim(observed)[1L])) {
          direct <- fit$main_effect$direct
          observed[t, , ] <- matrix(direct$c_hat[[m]][t, 1L], dim(observed)[2L], dim(observed)[3L]) +
            outer(direct$a_hat[[m]][t, ], rep(1, dim(observed)[3L])) +
            outer(rep(1, dim(observed)[2L]), direct$b_hat[[m]][t, ]) + fit$check_Y[[m]][t, , ]
        }
        expect_equal(observed, x$Xt[[m]], tolerance = 1e-12)
        expect_identical(dim(errors[[m]]), dim(x$Xt[[m]]))
        expect_equal(errors[[m]], x$Xt[[m]] - fitted(fit)[[m]], tolerance = 1e-12)
        expect_equal(errors[[m]], array(0, dim(errors[[m]])), tolerance = 1e-10)
      }
    }
  }
})

test_that("automatic and explicit ranks use the same final estimation path", {
  x <- .fit_fixture(rep(3L, 3L), rep(3L, 3L), c(1, 1, -2), c(1, 0, -1))
  rank <- select_MMEFM_rank(x$Xt, K0 = 8L, seed = 2026L)
  explicit <- est_MMEFM(x$Xt, rank, K0 = 8L, seed = 2026L)
  automatic <- est_MMEFM(x$Xt, K0 = 8L, seed = 2026L)
  expect_identical(explicit[-1L], automatic[-1L])
  expect_identical(explicit$rank$diagnostics, rank$diagnostics)
  reordered <- rank
  for (name in c("r2", "l2", "kr_m", "kc_m")) reordered[[name]] <- as.numeric(rank[[name]])[3:1]
  for (name in c("r2", "l2", "kr_m", "kc_m")) names(reordered[[name]]) <- names(x$Xt)[3:1]
  canonical <- est_MMEFM(x$Xt, reordered, K0 = 8L)
  expect_identical(canonical$rank, rank)
})

test_that("estimator dimension checks are method-specific", {
  Xt <- list(a = array(0, c(1L, 4L, 4L)), b = array(0, c(1L, 5L, 5L)))
  rank <- list(r1 = 1L, l1 = 1L, r2 = c(1L, 1L), l2 = c(1L, 1L), kr = 1L, kc = 1L, kr_m = c(1L, 1L), kc_m = c(1L, 1L))
  for (name in c("r1", "l1")) {
    bad <- rank
    bad[[name]] <- 2L
    expect_error(est_MMEFM(Xt, bad), "must not exceed T")
  }
  Xt <- lapply(Xt, function(a) array(0, c(2L, dim(a)[2:3])))
  rank$kr <- rank$kc <- 2L
  expect_error(est_MMEFM(Xt, rank, refit_method = "ALS"), "kr \\* kc <= T")
  expect_error(est_MMEFM(Xt, rank, alignment_method = "Procrustes"), "kr \\* kc <= T")
  expect_error(est_MMEFM(Xt, rank, alignment_method = "Procrustes", refit_method = "ALS"), "singular")
})

test_that("public estimator controls are validated", {
  x <- .fit_fixture()
  invalid <- list(K0 = 0, max_iter = 1.5, tol = Inf, max_iter_procrustes = 0, tol_procrustes = NA_real_,
                  max_iter_als = 0, tol_als = 0, lambda = -1, seed = -1, verbose = NA,
                  alignment_method = "other", refit_method = "other")
  for (name in names(invalid)) {
    args <- c(list(Xt = x$Xt, rank = x$rank), setNames(list(invalid[[name]]), name))
    expect_error(do.call(est_MMEFM, args), if (name %in% c("alignment_method", "refit_method")) "arg" else name)
  }
})

test_that("estimation preserves caller RNG on return and stochastic error", {
  x <- .fit_fixture()
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (existed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(23L)
  saved <- .Random.seed
  first <- est_MMEFM(x$Xt, x$rank, K0 = 8L)
  expect_identical(.Random.seed, saved)
  set.seed(71L)
  second <- est_MMEFM(x$Xt, x$rank, K0 = 8L)
  expect_identical(first[-1L], second[-1L])
  saved <- .Random.seed
  bad <- x$rank
  bad$kr <- bad$kc <- 2L
  expect_error(est_MMEFM(x$Xt, bad, K0 = 8L), "singular")
  expect_identical(.Random.seed, saved)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(est_MMEFM(x$Xt, x$rank, K0 = 8L))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("iteration limits warn once for each active stage", {
  x <- .fit_fixture()
  set.seed(91L)
  noisy <- lapply(x$Xt, function(a) a + array(rnorm(length(a), sd = 0.05), dim(a)))
  cases <- list(list(max_iter = 1L, tol = 1e-15), list(alignment_method = "Procrustes", max_iter_procrustes = 1L), list(refit_method = "ALS", max_iter_als = 1L))
  patterns <- c("Common-loading projection", "Procrustes.*a, b", "ALS refitting")
  for (i in seq_along(cases)) {
    messages <- character()
    withCallingHandlers(do.call(est_MMEFM, c(list(Xt = noisy, rank = x$rank, K0 = 8L), cases[[i]])), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_length(messages, 1L)
    expect_match(messages, patterns[i])
  }
  a <- .fit_fixture(rep(3L, 3L), rep(3L, 3L), c(1, 1, -2), c(1, 0, -1))
  messages <- character()
  withCallingHandlers(select_MMEFM_rank(a$Xt, K0 = 8L, max_iter = 1L, tol = 1e-20), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(messages, 1L)
  expect_match(messages, "Rank-selection projection")
})

test_that("print and summary show only compact structure", {
  x <- .fit_fixture()
  fit <- est_MMEFM(x$Xt, x$rank, K0 = 8L)
  result <- summary(fit)
  expect_s3_class(result, "summary.mmefm_fit")
  expect_identical(names(result), c("call", "dimensions", "rank", "alignment_method", "refit_method"))
  expect_identical(names(result$rank), names(x$rank))
  output <- capture.output(returned <- withVisible(print(fit)))
  expect_false(returned$visible)
  expect_identical(returned$value, fit)
  expect_lt(length(output), 12L)
  expect_match(paste(output, collapse = " "), "Global ranks.*r1=1.*VanLoan")
  expect_false(any(grepl("converg|iterations|eigenvalues|ratios|check_Y", output)))
  expect_identical(capture.output(print(result)), output)
})

test_that("structured main reconstruction differs from direct moments under ridge", {
  x <- .fit_fixture()
  fit <- est_MMEFM(x$Xt, x$rank, K0 = 8L, lambda = 0.1)
  direct_main <- Map(`-`, x$Xt, fit$check_Y)
  expect_gt(max(abs(unlist(Map(`-`, fitted(fit, "main"), direct_main)))), 1e-4)
  sum_components <- Map(`+`, Map(`+`, fitted(fit, "main"), fitted(fit, "global_common")), fitted(fit, "local_common"))
  expect_equal(fitted(fit), sum_components, tolerance = 1e-12)
})

test_that("fitted and residual methods retain a singleton time dimension", {
  x <- .fit_fixture()
  fit <- est_MMEFM(x$Xt, x$rank, K0 = 8L)
  fit$dimensions$T <- 1L
  for (name in c("mu_hat", "alpha_hat", "beta_hat")) fit$main_effect$global[[name]] <- fit$main_effect$global[[name]][1L, , drop = FALSE]
  for (name in c("mu_hat", "alpha_hat", "beta_hat")) {
    fit$main_effect$local[[name]] <- lapply(fit$main_effect$local[[name]], function(a) a[1L, , drop = FALSE])
  }
  for (name in c("c_hat", "a_hat", "b_hat")) fit$main_effect$direct[[name]] <- lapply(fit$main_effect$direct[[name]], function(a) a[1L, , drop = FALSE])
  fit$common_component$global$G_tilde <- fit$common_component$global$G_tilde[1L, , , drop = FALSE]
  fit$common_component$local$F_hat <- lapply(fit$common_component$local$F_hat, function(a) a[1L, , , drop = FALSE])
  fit$check_Y <- lapply(fit$check_Y, function(a) a[1L, , , drop = FALSE])
  for (component in c("all", "main", "global_main", "local_main", "global_common", "local_common")) {
    values <- fitted(fit, component)
    for (m in seq_along(values)) expect_identical(dim(values[[m]]), c(1L, dim(x$Xt[[m]])[2:3]))
  }
  for (m in seq_along(x$Xt)) expect_equal(residuals(fit)[[m]], array(0, c(1L, dim(x$Xt[[m]])[2:3])), tolerance = 1e-10)
})
