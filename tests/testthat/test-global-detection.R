test_that("RMS standardization uses the supplied residuals without floors", {
  Y <- list(a = array(seq_len(24L), c(2L, 3L, 4L)), b = array(seq_len(30L) / 2, c(2L, 5L, 3L)))
  normalized <- .standardize_detection_residuals(Y)
  expect_identical(names(normalized$check_Z), names(Y))
  for (m in seq_along(Y)) {
    expect_equal(unname(normalized$scale[m]^2), mean(Y[[m]]^2), tolerance = 1e-14)
    expect_equal(normalized$check_Z[[m]], Y[[m]] / normalized$scale[m], tolerance = 1e-14)
  }
  Y$a[] <- 0
  expect_error(.standardize_detection_residuals(Y), "scale.*group a")
  Y$a[] <- 1e200
  expect_error(.standardize_detection_residuals(Y), "scale.*group a")
})

test_that("S_G follows the entrywise manuscript sums and maxima", {
  Y <- list(a = array(sin(seq_len(48L)) + 1, c(4L, 3L, 4L)),
            b = array(cos(seq_len(60L)) + 2, c(4L, 5L, 3L)),
            c = array(sin(seq_len(64L) / 3) - 1, c(4L, 4L, 4L)))
  Y <- .standardize_detection_residuals(Y)$check_Z
  rank <- list(kr = 2L, kc = 2L, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L))
  Q <- J <- vector("list", 3L)
  for (m in seq_along(Y)) {
    Q[[m]] <- diag(dim(Y[[m]])[2L])[, 1:2, drop = FALSE]
    J[[m]] <- diag(dim(Y[[m]])[3L])[, 1:2, drop = FALSE]
  }
  result <- .global_detection_statistic(Y, rank, Q, J)
  expected <- numeric(3L)
  for (m in seq_along(Y)) {
    numerator <- 0
    for (n in seq_along(Y)) {
      if (n == m) next
      covariances <- numeric()
      for (cm in seq_len(rank$kr)) for (cn in seq_len(rank$kr)) for (dm in seq_len(rank$kc)) for (dn in seq_len(rank$kc)) {
        covariance <- 0
        for (t in seq_len(4L)) {
          covariance <- covariance + as.numeric(t(Q[[m]][, cm]) %*% Y[[m]][t, , ] %*% J[[m]][, dm]) *
            as.numeric(t(Q[[n]][, cn]) %*% Y[[n]][t, , ] %*% J[[n]][, dn]) / 4
        }
        covariances <- c(covariances, covariance^2)
      }
      numerator <- numerator + max(covariances)
    }
    row_cov <- matrix(0, dim(Y[[m]])[2L], dim(Y[[m]])[2L])
    col_cov <- matrix(0, dim(Y[[m]])[3L], dim(Y[[m]])[3L])
    for (t in seq_len(4L)) {
      row_cov <- row_cov + Y[[m]][t, , ] %*% t(Y[[m]][t, , ]) / 4
      col_cov <- col_cov + t(Y[[m]][t, , ]) %*% Y[[m]][t, , ] / 4
    }
    U <- eigen((row_cov + t(row_cov)) / 2, symmetric = TRUE)$vectors[, 1:3, drop = FALSE]
    V <- eigen((col_cov + t(col_cov)) / 2, symmetric = TRUE)$vectors[, 1:3, drop = FALSE]
    energies <- numeric()
    for (c in seq_len(3L)) for (d in seq_len(3L)) {
      energy <- 0
      for (t in seq_len(4L)) energy <- energy + as.numeric(t(U[, c]) %*% Y[[m]][t, , ] %*% V[, d])^2 / 4
      energies <- c(energies, energy)
    }
    expected[m] <- numerator / max(energies)^2
  }
  names(expected) <- names(Y)
  expect_equal(result$statistic, expected, tolerance = 1e-13)
  expect_equal(result$second_largest, unname(sort(expected, decreasing = TRUE)[2L]), tolerance = 1e-13)
  zero <- Y
  zero$a[] <- 0
  expect_error(.global_detection_statistic(zero, rank, Q, J), "denominator.*group a")
})

test_that("circular shifts keep the reference and whole panels intact", {
  Y <- list(a = array(seq_len(24L), c(4L, 2L, 3L)), b = array(seq_len(32L), c(4L, 4L, 2L)), c = array(seq_len(36L), c(4L, 3L, 3L)))
  set.seed(2026L)
  shifts <- c(sample.int(4L, 1L) - 1L, sample.int(4L, 1L) - 1L)
  set.seed(2026L)
  shifted <- .circular_shift_detection_groups(Y, 2L)
  expect_identical(shifted$b, Y$b)
  expect_identical(names(shifted), names(Y))
  for (i in seq_along(shifts)) {
    m <- c(1L, 3L)[i]
    index <- 1L + ((seq_len(4L) - 1L + shifts[i]) %% 4L)
    expect_identical(shifted[[m]], Y[[m]][index, , , drop = FALSE])
    expect_identical(dim(shifted[[m]]), dim(Y[[m]]))
  }
  singleton <- lapply(Y, function(a) a[1L, , , drop = FALSE])
  expect_identical(.circular_shift_detection_groups(singleton, 2L), singleton)
})

test_that("raw and fit inputs share standardized re-estimation and exact contracts", {
  x <- .fit_fixture()
  fit <- est_MMEFM(x$Xt, x$rank, K0 = 8L)
  expect_warning(raw <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 8L), NA)
  from_fit <- detect_MMEFM_global(fit, B = 3L, K0 = 8L)
  expect_identical(raw[-1L], from_fit[-1L])
  expect_s3_class(raw, "mmefm_global_detection")
  expect_identical(names(raw), c("call", "statistic", "second_largest", "cutoff", "exists", "initial_selected", "selected", "scale", "rank", "bootstrap"))
  expect_identical(names(raw$bootstrap), c("B", "alpha", "reference_group", "second_largest"))
  expect_s3_class(raw$rank, "mmefm_rank")
  expect_null(raw$rank$diagnostics)
  expect_length(raw$bootstrap$second_largest, 3L)
  for (name in c("statistic", "initial_selected", "selected", "scale")) expect_identical(names(raw[[name]]), names(x$Xt))
  expect_type(raw$selected, "logical")
  expect_identical(raw$cutoff, quantile(raw$bootstrap$second_largest, 1 - raw$bootstrap$alpha, type = 8L, names = FALSE))
  expect_identical(raw$exists, raw$second_largest > raw$cutoff)
  expect_identical(raw$initial_selected, raw$statistic > raw$cutoff)
  expected <- raw$initial_selected
  if (!raw$exists || sum(expected) < 2L) expected[] <- FALSE
  expect_identical(raw$selected, expected)
  changed <- fit
  for (name in c("Q_hat", "J_hat", "Q_tilde", "J_tilde")) changed$common_component$global[[name]] <- lapply(changed$common_component$global[[name]], function(a) a + 100)
  changed$dimensions$T <- 999L
  expect_identical(detect_MMEFM_global(changed, x$rank, B = 3L, K0 = 8L)[-1L], raw[-1L])
  alternative <- x$rank
  alternative$kr_m <- alternative$kc_m <- rep(2L, 3L)
  expect_identical(detect_MMEFM_global(fit, alternative, B = 2L, K0 = 8L)$rank$kr_m, c(a = 2L, b = 2L, c = 2L))
})

test_that("automatic ranks and supplied ranks give identical detection", {
  x <- .fit_fixture(rep(3L, 3L), rep(3L, 3L), c(1, 1, -2), c(1, -1, 0))
  rank <- select_MMEFM_rank(x$Xt, K0 = 8L)
  explicit <- detect_MMEFM_global(x$Xt, rank, B = 3L, K0 = 8L)
  automatic <- detect_MMEFM_global(x$Xt, B = 3L, K0 = 8L)
  expect_identical(explicit[-1L], automatic[-1L])
  expect_identical(explicit$rank$diagnostics, rank$diagnostics)
})

test_that("detection restores caller RNG on success and stochastic error", {
  x <- .fit_fixture()
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (existed) old_seed <- .Random.seed
  on.exit({
    if (existed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(71L)
  saved <- .Random.seed
  first <- detect_MMEFM_global(x$Xt, x$rank, B = 2L, K0 = 8L)
  expect_identical(.Random.seed, saved)
  set.seed(92L)
  second <- detect_MMEFM_global(x$Xt, x$rank, B = 2L, K0 = 8L)
  expect_identical(first[-1L], second[-1L])
  saved <- .Random.seed
  rank <- .validate_rank(x$rank, .validate_Xt(x$Xt))
  Y <- .standardize_detection_residuals(.main_effect_moments(x$Xt)$check_Y)$check_Z
  Y$a[] <- 0
  expect_error(.with_preserved_seed(2026L, .bootstrap_global_detection(Y, rank, 2L, 1L, 8L, 20L, 1e-4, FALSE)), "Bootstrap replication 1 failed:.*denominator.*a")
  expect_identical(.Random.seed, saved)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(detect_MMEFM_global(x$Xt, x$rank, B = 2L, K0 = 8L))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("public control and scale errors are clear", {
  x <- .fit_fixture()
  invalid <- list(B = 0L, alpha = 1, reference_group = 4L, K0 = 1.5, max_iter = 0L, tol = Inf, seed = -1, verbose = NA)
  for (name in names(invalid)) {
    args <- list(Xt = x$Xt, rank = x$rank, B = 3L, K0 = 8L)
    args[[name]] <- invalid[[name]]
    expect_error(do.call(detect_MMEFM_global, args), name)
  }
  expect_error(detect_MMEFM_global(x$Xt, list(r1 = 1L), B = 3L), "complete named list")
  expect_error(detect_MMEFM_global(list(matrix(1, 3L, 3L)), x$rank, B = 3L), "at least two")
  zero <- lapply(x$Xt, function(a) array(0, dim(a)))
  expect_error(detect_MMEFM_global(zero, x$rank, B = 3L), "scale.*group a")
})

test_that("projection limits produce one observed and one aggregate warning", {
  x <- .fit_fixture()
  set.seed(91L)
  Y <- lapply(x$Xt, function(a) a + array(rnorm(length(a), sd = 0.05), dim(a)))
  messages <- character()
  withCallingHandlers(detect_MMEFM_global(Y, x$rank, B = 3L, K0 = 8L, max_iter = 1L, tol = 1e-15), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(messages, 2L)
  expect_match(messages[1L], "Observed detection projection")
  expect_match(messages[2L], "[1-3] of 3 bootstrap projections")
})

test_that("detection printing is compact and invisible", {
  x <- .fit_fixture()
  result <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 8L)
  output <- capture.output(returned <- withVisible(print(result)))
  expect_false(returned$visible)
  expect_identical(returned$value, result)
  expect_length(output, 4L)
  expect_match(paste(output, collapse = " "), "cutoff.*detected.*Selected groups")
  expect_false(any(grepl("bootstrap|eigenvalues|ratios", output)))
  result$selected[] <- FALSE
  expect_match(paste(capture.output(print(result)), collapse = " "), "Selected groups: none")
  result$selected[] <- TRUE
  expect_match(paste(capture.output(print(result)), collapse = " "), "Selected groups: a, b, c")
  names(result$selected) <- NULL
  expect_match(paste(capture.output(print(result)), collapse = " "), "Selected groups: 1, 2, 3")
})

test_that("strict cutoff ties give an empty screen with singleton time and unnamed groups", {
  Y <- rep(list(array(tcrossprod(c(1, -1, 0)), c(1L, 3L, 3L))), 3L)
  rank <- list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3L), l2 = rep(1L, 3L), kr = 1L, kc = 1L, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L))
  result <- detect_MMEFM_global(Y, rank, B = 3L, K0 = 2L)
  expect_equal(result$second_largest, result$cutoff, tolerance = 1e-14)
  expect_false(result$exists)
  expect_identical(result$initial_selected, rep(FALSE, 3L))
  expect_identical(result$selected, rep(FALSE, 3L))
  expect_null(names(result$scale))
  expect_null(names(result$statistic))
})

test_that("bootstrap draws follow one stream and re-estimate at fixed ranks", {
  x <- .fit_fixture()
  rank <- .validate_rank(x$rank, .validate_Xt(x$Xt))
  Z <- .standardize_detection_residuals(.main_effect_moments(x$Xt)$check_Y)$check_Z
  expected <- .with_preserved_seed(2026L, {
    invisible(.initial_common_loadings(Z, rank$kr, rank$kc, K0 = 8L))
    draws <- numeric(3L)
    for (b in seq_len(3L)) {
      shifted <- Z
      for (m in 2:3) {
        shift <- sample.int(16L, 1L) - 1L
        index <- 1L + ((seq_len(16L) - 1L + shift) %% 16L)
        shifted[[m]] <- Z[[m]][index, , , drop = FALSE]
      }
      loadings <- .initial_common_loadings(shifted, rank$kr, rank$kc, K0 = 8L)
      draws[b] <- .global_detection_statistic(shifted, rank, loadings$Q_hat, loadings$J_hat)$second_largest
    }
    draws
  })
  result <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 8L)
  expect_identical(result$bootstrap$second_largest, expected)
})

test_that("precomputed shifts and half-panels reproduce the sampled path without RNG", {
  x <- .fit_fixture()
  Z <- .standardize_detection_residuals(.main_effect_moments(x$Xt)$check_Y)$check_Z
  set.seed(2026L)
  shifts <- integer(length(Z))
  for (m in seq_along(Z)) if (m != 2L) shifts[m] <- sample.int(16L, 1L) - 1L
  indices <- vector("list", length(Z))
  for (m in seq_along(Z)) {
    q <- dim(Z[[m]])[3L]
    indices[[m]] <- vector("list", 4L)
    for (k in seq_len(4L)) indices[[m]][[k]] <- sample(q, floor(q / 2))
  }
  saved <- .Random.seed
  shifted <- .circular_shift_detection_groups(Z, 2L, shifts)
  expect_identical(.Random.seed, saved)
  expect_identical(shifted$b, Z$b)
  for (m in seq_along(Z)) {
    index <- 1L + ((seq_len(16L) - 1L + shifts[m]) %% 16L)
    expect_identical(shifted[[m]], Z[[m]][index, , , drop = FALSE])
  }
  expect_identical(.circular_shift_detection_groups(Z, 2L, integer(3L)), Z)
  fixed <- .initial_common_loadings(shifted, 1L, 1L, K0 = 4L, initial_candidate_indices = indices)
  expect_identical(.Random.seed, saved)
  inputs <- list(check_Z = Z, rank = x$rank, reference_group = 2L, K0 = 4L, max_iter = 20L, tol = 1e-4)
  draw <- .bootstrap_global_detection_one(list(b = 1L, shifts = shifts, initial_candidate_indices = indices), inputs)
  expect_identical(.Random.seed, saved)
  set.seed(2026L)
  sampled_shifted <- .circular_shift_detection_groups(Z, 2L)
  sampled <- .initial_common_loadings(sampled_shifted, 1L, 1L, K0 = 4L)
  expect_identical(sampled_shifted, shifted)
  expect_identical(sampled, fixed)
  expect_identical(draw$second_largest, .global_detection_statistic(shifted, x$rank, sampled$Q_hat, sampled$J_hat)$second_largest)
  expect_identical(draw$converged, sampled$converged)
})

test_that("parallel controls are validated even for serial execution", {
  x <- .fit_fixture()
  for (value in list(NA, 1, c(TRUE, FALSE))) {
    expect_error(detect_MMEFM_global(x$Xt, x$rank, B = 2L, parallel = value), "parallel")
  }
  for (value in list(0, -1, 1.5, NA, Inf, .Machine$integer.max + 1)) {
    expect_error(detect_MMEFM_global(x$Xt, x$rank, B = 2L, parallel = FALSE, num.cores = value), "num.cores")
  }
})

test_that("PSOCK bootstrap exactly matches serial results and preserves caller RNG", {
  x <- .fit_fixture()
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (existed) old_seed <- .Random.seed
  on.exit({
    if (existed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(71L)
  saved <- .Random.seed
  kind <- RNGkind()
  serial <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 4L, reference_group = 2L)
  parallel_fit <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 4L, reference_group = 2L, parallel = TRUE, num.cores = 2L)
  expect_identical(.Random.seed, saved)
  expect_identical(RNGkind(), kind)
  expect_identical(parallel_fit[-1L], serial[-1L])
  expect_identical(class(parallel_fit), class(serial))
  expect_identical(names(parallel_fit), names(serial))
  expect_identical(names(parallel_fit$bootstrap), names(serial$bootstrap))
  rm(".Random.seed", envir = .GlobalEnv)
  again <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 4L, reference_group = 2L, parallel = TRUE, num.cores = 2L)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_identical(again[-1L], serial[-1L])
})

test_that("one worker and one draw retain serial execution", {
  x <- .fit_fixture()
  serial <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 4L)
  one_worker <- detect_MMEFM_global(x$Xt, x$rank, B = 3L, K0 = 4L, parallel = TRUE, num.cores = 1)
  expect_identical(one_worker[-1L], serial[-1L])
  serial <- detect_MMEFM_global(x$Xt, x$rank, B = 1L, K0 = 4L)
  one_draw <- detect_MMEFM_global(x$Xt, x$rank, B = 1L, K0 = 4L, parallel = TRUE, num.cores = 2L)
  expect_identical(one_draw[-1L], serial[-1L])
})

test_that("worker numerical failures identify the draw and clean up sockets and RNG", {
  x <- .fit_fixture()
  rank <- .validate_rank(x$rank, .validate_Xt(x$Xt))
  Z <- .standardize_detection_residuals(.main_effect_moments(x$Xt)$check_Y)$check_Z
  # The same degenerate denominator fixture used by the serial failure regression.
  Z$a[] <- 0
  set.seed(71L)
  saved <- .Random.seed
  connections <- showConnections(all = TRUE)
  expect_error(.with_preserved_seed(2026L,
    .bootstrap_global_detection(Z, rank, 2L, 1L, 4L, 20L, 1e-4, FALSE, TRUE, 2L)),
    "Bootstrap replication 1 failed:.*denominator.*a")
  expect_identical(.Random.seed, saved)
  expect_identical(showConnections(all = TRUE), connections)
})

test_that("parallel projection limits produce only observed and aggregate warnings", {
  x <- .fit_fixture()
  set.seed(91L)
  Y <- lapply(x$Xt, function(a) a + array(rnorm(length(a), sd = 0.05), dim(a)))
  messages <- character()
  withCallingHandlers(detect_MMEFM_global(Y, x$rank, B = 2L, K0 = 4L, max_iter = 1L,
                                          tol = 1e-15, parallel = TRUE, num.cores = 2L), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(messages, 2L)
  expect_match(messages[1L], "Observed detection projection")
  expect_match(messages[2L], "[1-2] of 2 bootstrap projections")
})
