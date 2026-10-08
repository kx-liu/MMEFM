test_that("piecewise series retain matrix dimensions and canonical segment metadata", {
  expect_identical(dim(gen_piecewise_arfima(3L)), c(3L, 1L))
  shared <- gen_piecewise_arfima(7L, 2L, c(-0.8, 0.2, 0), burn = 2L,
                               filter_length = 4L, return_filters = TRUE)
  expect_identical(dim(shared$series), c(7L, 2L))
  expect_true(all(is.finite(shared$series)))
  expect_identical(names(shared), c("series", "phi", "filter_length", "short_filter", "series_specs"))
  expect_identical(shared$series_specs[[1L]], shared$series_specs[[2L]])
  spec <- shared$series_specs[[1L]]
  expect_identical(names(spec), c("d", "seg_lens", "seg_start", "seg_end", "fractional_filters"))
  expect_identical(spec$d, c(-0.8, 0.2, 0))
  expect_identical(spec$seg_lens, c(2L, 2L, 3L))
  expect_identical(spec$seg_start, c(1L, 3L, 5L))
  expect_identical(spec$seg_end, c(2L, 4L, 7L))
  expect_identical(shared$short_filter, 0.5^(0:4))
  for (j in seq_along(spec$d)) {
    expected <- numeric(5L)
    expected[1L] <- 1
    for (k in seq_len(4L)) {
      expected[k + 1L] <- expected[k] * (k - 1 + spec$d[j]) / k
    }
    expect_identical(spec$fractional_filters[[j]], expected)
  }
  explicit <- list(d = c(0.1, -0.2), seg_lens = c(2L, 5L))
  x <- gen_piecewise_arfima(7L, 2L, explicit, burn = 2L, filter_length = 4L, return_filters = TRUE)
  expect_identical(x$series_specs[[1L]]$seg_lens, c(2L, 5L))
  individual <- gen_piecewise_arfima(7L, 2L, list(c(0, 0.2), explicit),
                                   burn = 2L, filter_length = 4L, return_filters = TRUE)
  expect_identical(dim(individual$series), c(7L, 2L))
  expect_identical(individual$series_specs[[1L]]$seg_lens, c(3L, 4L))
  expect_identical(individual$series_specs[[2L]]$seg_lens, c(2L, 5L))
  expect_identical(individual, gen_piecewise_arfima(7L, 2L, list(c(0, 0.2), explicit),
                                                  burn = 2L, filter_length = 4L, return_filters = TRUE))
})

test_that("zero-memory AR(p) recursion and raw innovations agree with explicit calculations", {
  for (family in c("Gaussian", "Student-t")) {
    phi <- c(0.4, -0.2)
    innovations <- .with_preserved_seed(2026L, {
      if (family == "Gaussian") matrix(stats::rnorm(16L), 8L, 2L) else matrix(stats::rt(16L, df = 1), 8L, 2L)
    })
    expected <- innovations
    for (t in seq.int(2L, 8L)) {
      for (lag in seq_len(min(2L, t - 1L))) {
        expected[t, ] <- expected[t, ] + phi[lag] * expected[t - lag, ]
      }
    }
    x <- gen_piecewise_arfima(8L, 2L, phi = phi, innovation = family, df = 1,
                             burn = 0L, filter_length = 0L)
    expect_identical(x, expected)
    raw <- gen_piecewise_arfima(8L, 2L, phi = 0, innovation = family, df = 1,
                               burn = 0L, filter_length = 0L)
    expect_identical(raw, innovations)
  }
})

test_that("nonzero-memory AR(p) filters have the specified impulse response and normalization", {
  phi <- c(1.2, -0.5)
  x <- gen_piecewise_arfima(6L, 2L, c(0.1, -0.2), phi = phi,
                           burn = 2L, filter_length = 4L, return_filters = TRUE)
  expect_equal(x$short_filter, c(1, 1.2, 0.94, 0.528, 0.1636), tolerance = 1e-14)
  innovations <- .with_preserved_seed(2026L, matrix(stats::rnorm(32L), 16L, 2L))
  short <- as.matrix(stats::filter(innovations, x$short_filter, sides = 1L))
  expected <- matrix(0, 6L, 2L)
  spec <- x$series_specs[[1L]]
  for (j in seq_along(spec$d)) {
    block <- short[seq.int(10L + spec$seg_start[j] - 4L, 10L + spec$seg_end[j]), , drop = FALSE]
    filtered <- as.matrix(stats::filter(block, spec$fractional_filters[[j]], sides = 1L))
    expected[seq.int(spec$seg_start[j], spec$seg_end[j]), ] <-
      filtered[seq.int(nrow(filtered) - spec$seg_lens[j] + 1L, nrow(filtered)), , drop = FALSE]
  }
  expect_identical(x$series, expected)
  expect_true(all(is.finite(x$series)))
  expect_identical(.simulation_ar_filter(c(0, 0), 4L), c(1, 0, 0, 0, 0))
  expect_identical(.simulation_ar_filter(c(0.5, 0), 4L), 0.5^(0:4))
  expect_identical(dim(gen_piecewise_arfima(1L, phi = c(0, 0), burn = 0L, filter_length = 0L)), c(1L, 1L))
})

test_that("invalid AR coefficients and memory specifications error without repair", {
  for (phi in list(numeric(), "0.5", matrix(0.5), c(0.2, NA), Inf, 1, -1,
                   c(1, 0), c(0.5, 0.5), c(0, -1), c(1.4, 0.1))) {
    expect_error(gen_piecewise_arfima(6L, phi = phi, burn = 0L, filter_length = 2L), "phi")
  }
  for (spec in list(0.5, Inf, numeric(), rep(0, 7L),
                    list(d = c(0, 0.1), seg_lens = c(2L, 3L)),
                    list(d = c(0, 0.1), seg_lens = c(0L, 6L)),
                    list(d = c(0, 0.1), seg_lens = c(2.5, 3.5)),
                    list(d = c(0, 0.1), seg_lens = c(Inf, 6L)), list(0, 0))) {
    expect_error(gen_piecewise_arfima(6L, d_vec = spec, burn = 0L, filter_length = 2L), "d_vec")
  }
  for (name in c("T", "n", "burn", "filter_length", "seed")) {
    for (value in list(-1, 1.5, NA_real_, Inf, numeric(), matrix(2))) {
      args <- list(T = 6L, burn = 0L, filter_length = 2L)
      args[[name]] <- value
      expect_error(do.call(gen_piecewise_arfima, args), name)
    }
  }
  expect_error(gen_piecewise_arfima(0L), "T")
  expect_error(gen_piecewise_arfima(6L, n = 0L), "n")
  expect_error(gen_piecewise_arfima(6L, df = 0), "df")
  expect_error(gen_piecewise_arfima(6L, return_filters = NA), "return_filters")
})

test_that("scoped RNG and filter metadata preserve exactly the same series", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- .Random.seed
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(71L)
  saved <- .Random.seed
  kind <- RNGkind()
  for (memory in list(0, c(0.1, 0.2), list(0, c(0.1, 0.2)))) {
    plain <- gen_piecewise_arfima(6L, 2L, memory, phi = c(0.4, -0.1), burn = 2L, filter_length = 4L)
    detailed <- gen_piecewise_arfima(6L, 2L, memory, phi = c(0.4, -0.1), burn = 2L,
                                   filter_length = 4L, return_filters = TRUE)
    expect_identical(plain, detailed$series)
    expect_identical(.Random.seed, saved)
    expect_identical(RNGkind(), kind)
  }
  expect_error(gen_piecewise_arfima(6L, d_vec = -1e308, burn = 2L, filter_length = 4L))
  expect_identical(.Random.seed, saved)
  expect_identical(RNGkind(), kind)
  rm(".Random.seed", envir = .GlobalEnv)
  gen_piecewise_arfima(2L, burn = 0L, filter_length = 0L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_error(gen_piecewise_arfima(6L, phi = 1), "phi")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
