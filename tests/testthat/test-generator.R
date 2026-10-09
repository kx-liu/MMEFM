.generator_args <- function() {
  list(T = 8L, p = c(a = 4L, b = 5L, c = 6L), q = c(c = 4L, b = 6L, a = 5L),
       rank = list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3L), l2 = rep(1L, 3L), kr = 1L, kc = 1L, kr_m = rep(1L, 3L), kc_m = rep(1L, 3L)),
       burn = 4L, filter_length = 8L)
}

test_that("simulation has canonical dimensions, truth names, and centered loadings", {
  args <- .generator_args()
  args$rank <- .new_mmefm_rank(args$rank, diagnostics = list(unused = TRUE))
  x <- do.call(gen_MMEFM, args)
  expect_identical(names(x), c("Xt", "dimensions", "rank", "main_effect", "common_component", "error", "memory", "innovation", "components"))
  expect_null(attr(x, "class"))
  expect_identical(x$dimensions, list(M = 3L, T = 8L, p = c(a = 4L, b = 5L, c = 6L), q = c(a = 5L, b = 6L, c = 4L), group_names = c("a", "b", "c")))
  expect_s3_class(x$rank, "mmefm_rank")
  expect_null(x$rank$diagnostics)
  expect_null(x$components)
  expect_identical(x$innovation, list(family = "Gaussian", df = 6))
  expect_identical(names(x$main_effect), c("direct", "global", "local"))
  expect_identical(names(x$common_component), c("global", "local"))
  expect_identical(names(x$error), c("global", "local", "sigma"))
  lists <- c(x$main_effect$direct, x$main_effect$global[c("A1", "B1")], x$main_effect$local,
             x$common_component$global[c("Q", "J")], x$common_component$local,
             x$error$global[c("Q", "J")], x$error$local, list(x$error$sigma))
  for (values in lists) expect_identical(names(values), c("a", "b", "c"))
  for (m in seq_len(3L)) {
    expect_identical(dim(x$Xt[[m]]), c(8L, unname(x$dimensions$p[m]), unname(x$dimensions$q[m])))
    expect_true(all(is.finite(x$Xt[[m]])))
    for (values in c(x$main_effect$global[c("A1", "B1")], x$main_effect$local[c("A2", "B2")], x$common_component$global[c("Q", "J")], x$common_component$local[c("R", "C")])) {
      expect_equal(colMeans(values[[m]]), rep(0, ncol(values[[m]])), tolerance = 1e-14)
      expect_identical(ncol(values[[m]]), 1L)
    }
    expect_identical(dim(x$common_component$local$F[[m]]), c(8L, 1L, 1L))
    expect_identical(dim(x$error$local$F[[m]]), c(8L, 2L, 2L))
    expect_identical(dim(x$error$global$Q[[m]]), c(unname(x$dimensions$p[m]), 2L))
    expect_identical(dim(x$error$global$J[[m]]), c(unname(x$dimensions$q[m]), 2L))
  }
  expect_identical(dim(x$common_component$global$G), c(8L, 1L, 1L))
  expect_identical(dim(x$error$global$G), c(8L, 2L, 2L))
  args$p <- unname(args$p)
  args$q <- unname(args$q)
  expect_null(names(do.call(gen_MMEFM, args)$Xt))
})

test_that("main score projection and grand means enforce exact sample IC2", {
  x <- do.call(gen_MMEFM, .generator_args())
  global <- x$main_effect$global
  local <- x$main_effect$local
  expect_equal(Reduce(`+`, local$mu), matrix(0, 8L, 1L), tolerance = 1e-13)
  for (m in seq_len(3L)) {
    expect_equal(crossprod(global$alpha, local$alpha[[m]]), matrix(0, 1L, 1L), tolerance = 1e-12)
    expect_equal(crossprod(global$beta, local$beta[[m]]), matrix(0, 1L, 1L), tolerance = 1e-12)
    expect_equal(x$main_effect$direct$c[[m]], global$mu + local$mu[[m]])
    expect_equal(x$main_effect$direct$a[[m]], global$alpha %*% t(global$A1[[m]]) + local$alpha[[m]] %*% t(local$A2[[m]]))
    expect_equal(x$main_effect$direct$b[[m]], global$beta %*% t(global$B1[[m]]) + local$beta[[m]] %*% t(local$B2[[m]]))
  }
})

test_that("stored components recover Xt and structured errors independently", {
  args <- .generator_args()
  plain <- do.call(gen_MMEFM, args)
  args$store_components <- TRUE
  x <- do.call(gen_MMEFM, args)
  expect_identical(x[names(plain)[-length(plain)]], plain[-length(plain)])
  expect_identical(names(x$components), c("global_main", "local_main", "global_common", "local_common", "error_global", "error_local", "error_idiosyncratic"))
  for (values in x$components) expect_identical(names(values), names(x$Xt))
  for (m in seq_len(3L)) {
    terms <- lapply(x$components, function(values) values[[m]])
    expect_equal(Reduce(`+`, terms), x$Xt[[m]], tolerance = 1e-13)
    for (t in seq_len(8L)) {
      expect_equal(x$components$error_global[[m]][t, , ], x$error$global$Q[[m]] %*% matrix(x$error$global$G[t, , ], 2L, 2L) %*% t(x$error$global$J[[m]]), tolerance = 1e-13)
      expect_equal(x$components$error_local[[m]][t, , ], x$error$local$R[[m]] %*% matrix(x$error$local$F[[m]][t, , ], 2L, 2L) %*% t(x$error$local$C[[m]]), tolerance = 1e-13)
    }
    # The API retains sigma rather than duplicating the idiosyncratic innovations.
    other_terms <- terms[-length(terms)]
    epsilon <- sweep(x$Xt[[m]] - Reduce(`+`, other_terms), c(2L, 3L), x$error$sigma[[m]], "/")
    expect_equal(sweep(epsilon, c(2L, 3L), x$error$sigma[[m]], "*"), x$components$error_idiosyncratic[[m]], tolerance = 1e-12)
  }
  args$error_rank <- list(ker = 0L, kec = 2L, ker_m = c(1L, 0L, 0L), kec_m = c(0L, 0L, 2L))
  zero <- do.call(gen_MMEFM, args)
  expect_identical(dim(zero$error$global$G), c(8L, 0L, 2L))
  expect_identical(dim(zero$error$global$Q$a), c(4L, 0L))
  for (m in seq_len(3L)) {
    expect_true(all(zero$components$error_global[[m]] == 0))
    expect_true(all(zero$components$error_local[[m]] == 0))
    expect_identical(dim(zero$error$local$F[[m]]), c(8L, args$error_rank$ker_m[m], args$error_rank$kec_m[m]))
  }
})

test_that("ownership mutation preserves the full draw order", {
  args <- .generator_args()
  args$store_components <- TRUE
  all <- do.call(gen_MMEFM, args)
  args$active_global <- c(3L, 1L)
  pair <- do.call(gen_MMEFM, args)
  args$active_global <- integer(0L)
  null <- do.call(gen_MMEFM, args)
  for (x in list(pair, null)) {
    for (name in c("rank", "main_effect", "error", "memory", "innovation")) expect_identical(x[[name]], all[[name]])
    expect_identical(x$common_component$local, all$common_component$local)
    expect_identical(x$common_component$global$G, all$common_component$global$G)
    expect_identical(x$components[-3L], all$components[-3L])
    for (m in seq_len(3L)) {
      if (x$common_component$global$active[m]) {
        expect_identical(x$common_component$global$Q[[m]], all$common_component$global$Q[[m]])
        expect_identical(x$common_component$global$J[[m]], all$common_component$global$J[[m]])
      } else {
        expect_true(all(x$common_component$global$Q[[m]] == 0))
        expect_true(all(x$common_component$global$J[[m]] == 0))
        expect_true(all(x$components$global_common[[m]] == 0))
      }
    }
  }
  expect_identical(pair$common_component$global$active, c(a = TRUE, b = FALSE, c = TRUE))
  args$active_global <- c(TRUE, FALSE, TRUE)
  expect_identical(do.call(gen_MMEFM, args), pair)
  args$active_global <- 1L
  expect_error(do.call(gen_MMEFM, args), "singleton")
})

test_that("strength expansion and named memory specifications are canonical", {
  args <- .generator_args()
  args$rank$kr <- 2L
  args$strength <- list(common_global = list(list(0.8, 1), list(c(0.7, 0.9), 1), list(0.8, 1)))
  scalar <- do.call(gen_MMEFM, args)
  args$strength$common_global <- list(list(c(0.8, 0.8), 1), list(c(0.7, 0.9), 1), list(c(0.8, 0.8), 1))
  expect_identical(do.call(gen_MMEFM, args), scalar)
  args$memory <- list(mu_global = c(0.1, -0.1, 0), mu_local = c(0, 0),
                      alpha = list(global = list(d = 0.2, seg_lens = 8L), local = list(c = 0, b = c(0, 0.1), a = -0.1)),
                      beta = list(global = 0, local = c(0.1, 0)), G = 0,
                      F = list(c = c(0, -0.1), b = list(d = 0.1, seg_lens = 8L), a = 0))
  x <- do.call(gen_MMEFM, args)
  expect_identical(x$memory$mu_global, list(d = c(0.1, -0.1, 0), seg_lens = c(2L, 2L, 4L)))
  expect_identical(x$memory$alpha$local$a, list(d = -0.1, seg_lens = 8L))
  expect_identical(names(x$memory$F), c("a", "b", "c"))
  expect_identical(dim(x$common_component$global$Q$a), c(4L, 2L))
  for (spec in list(c(0, 0.5), list(d = 0, seg_lens = 7L), list(d = c(0, 0), seg_lens = c(3.5, 4.5)))) {
    args$memory$G <- spec
    expect_error(do.call(gen_MMEFM, args), "memory\\$G")
  }
})

test_that("temporal equality is allowed and infeasible rank sums error", {
  args <- .generator_args()
  args$T <- 2L
  expect_warning(do.call(gen_MMEFM, args), NA)
  args$T <- 1L
  expect_error(do.call(gen_MMEFM, args), "row.*group.*a")
  args$T <- 2L
  args$rank$l1 <- 2L
  expect_error(do.call(gen_MMEFM, args), "column.*group.*a")
  args$T <- 8L
  args$rank$r1 <- 3L
  args$rank$r2 <- rep(2L, 3L)
  expect_error(do.call(gen_MMEFM, args), "r1.*r2.*row")
})

test_that("innovations are raw Gaussian and raw Student-t draws", {
  set.seed(31L)
  expected <- stats::rt(12L, df = 2)
  set.seed(31L)
  expect_identical(.simulation_innovations(12L, "Student-t", 2), expected)
  set.seed(31L)
  expected <- stats::rnorm(12L)
  set.seed(31L)
  expect_identical(.simulation_innovations(12L, "Gaussian", 6), expected)
  args <- .generator_args()
  args$innovation <- "Student-t"
  args$df <- 2
  expect_true(all(is.finite(unlist(do.call(gen_MMEFM, args)$Xt))))
})

test_that("generator preserves RNG after success and stochastic error", {
  args <- .generator_args()
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (existed) old_seed <- .Random.seed
  on.exit({
    if (existed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(42L)
  saved <- .Random.seed
  first <- do.call(gen_MMEFM, args)
  expect_identical(.Random.seed, saved)
  set.seed(91L)
  saved <- .Random.seed
  expect_identical(do.call(gen_MMEFM, args), first)
  expect_identical(.Random.seed, saved)
  # An extreme finite negative d overflows its fractional filter after draws begin.
  args$memory <- setNames(rep(list(0), 6L), c("mu_global", "mu_local", "alpha", "beta", "G", "F"))
  args$memory$G <- -1e200
  expect_error(do.call(gen_MMEFM, args), "Nonfinite generated")
  expect_identical(.Random.seed, saved)
  args$memory <- NULL
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(do.call(gen_MMEFM, args))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("generator controls reject invalid dimensions and specifications", {
  args <- .generator_args()
  invalid <- list(T = 1.5, p = c(1, 4, 5), q = c(4, 5), strength = list(main_global = c(0.5, 1)), phi = 1, error_phi = Inf,
                  error_rank = list(ker = -1),
                  error_loading_zero_prob = 1, innovation = "other", df = 0, burn = -1,
                  filter_length = 1.5, active_global = c(1, 1), seed = NA_real_, store_components = NA)
  for (name in names(invalid)) {
    bad <- args
    bad[[name]] <- invalid[[name]]
    expect_error(do.call(gen_MMEFM, bad), if (name == "innovation") "arg" else if (name == "q") "same length" else name)
  }
  args$q <- unname(args$q)
  expect_error(do.call(gen_MMEFM, args), "both be named")
  args <- .generator_args()
  names(args$q) <- c("a", "a", "b")
  expect_error(do.call(gen_MMEFM, args), "unique")
  args <- .generator_args()
  args$memory <- list(G = 0)
  expect_error(do.call(gen_MMEFM, args), "exactly")
  args <- .generator_args()
  args$rank <- list(kr = 1L)
  expect_error(do.call(gen_MMEFM, args), "complete named list")
})

test_that("rectangular multi-rank cores and score matrices retain their shapes", {
  args <- .generator_args()
  args$rank$r1 <- args$rank$l1 <- args$rank$kr <- 2L
  args$rank$kc <- 3L
  args$q["c"] <- 5L
  args$rank$kr_m <- c(1L, 2L, 1L)
  args$store_components <- TRUE
  x <- do.call(gen_MMEFM, args)
  expect_identical(dim(x$main_effect$global$alpha), c(8L, 2L))
  expect_identical(dim(x$main_effect$global$beta), c(8L, 2L))
  expect_identical(dim(x$common_component$global$G), c(8L, 2L, 3L))
  expect_identical(dim(x$common_component$local$F$b), c(8L, 2L, 1L))
  for (m in seq_len(3L)) {
    expect_equal(Reduce(`+`, lapply(x$components, function(a) a[[m]])), x$Xt[[m]], tolerance = 1e-13)
    expect_equal(crossprod(x$main_effect$global$alpha, x$main_effect$local$alpha[[m]]), matrix(0, 2L, 1L), tolerance = 1e-12)
  }
})


test_that("maximum centered spatial ranks and the temporal boundary are usable", {
  p <- c(a = 4L, b = 5L, c = 6L)
  q <- c(a = 5L, b = 6L, c = 4L)
  rank <- list(r1 = 2L, l1 = 2L, r2 = c(a = 1L, b = 2L, c = 3L), l2 = c(a = 2L, b = 3L, c = 1L),
               kr = 2L, kc = 2L, kr_m = c(a = 1L, b = 2L, c = 3L), kc_m = c(a = 2L, b = 3L, c = 1L))
  x <- gen_MMEFM(5L, p, q, rank, burn = 4L, filter_length = 8L, seed = 2026L, store_components = TRUE)
  expect_identical(.validate_rank(x$rank, x$dimensions), rank)
  expect_identical(rank$r1 + rank$r2, p - 1L)
  expect_identical(rank$kr + rank$kr_m, p - 1L)
  expect_identical(rank$l1 + rank$l2, q - 1L)
  expect_identical(rank$kc + rank$kc_m, q - 1L)
  global <- x$main_effect$global
  local <- x$main_effect$local
  for (m in seq_along(p)) {
    expect_identical(dim(x$Xt[[m]]), c(5L, unname(p[m]), unname(q[m])))
    expect_true(all(is.finite(x$Xt[[m]])))
    expect_equal(Reduce(`+`, lapply(x$components, function(values) values[[m]])), x$Xt[[m]], tolerance = 1e-13)
    loadings <- list(cbind(global$A1[[m]], local$A2[[m]]), cbind(global$B1[[m]], local$B2[[m]]),
                     cbind(x$common_component$global$Q[[m]], x$common_component$local$R[[m]]),
                     cbind(x$common_component$global$J[[m]], x$common_component$local$C[[m]]))
    dimensions <- c(p[m], q[m], p[m], q[m])
    for (i in seq_along(loadings)) {
      expect_equal(colMeans(loadings[[i]]), rep(0, ncol(loadings[[i]])), tolerance = 1e-14)
      expect_identical(qr(loadings[[i]])$rank, unname(dimensions[i] - 1L))
    }
    expect_equal(crossprod(global$alpha, local$alpha[[m]]), matrix(0, rank$r1, rank$r2[m]), tolerance = 1e-12)
    expect_equal(crossprod(global$beta, local$beta[[m]]), matrix(0, rank$l1, rank$l2[m]), tolerance = 1e-12)
    expect_identical(qr(cbind(global$alpha, local$alpha[[m]]))$rank, unname(rank$r1 + rank$r2[m]))
    expect_identical(qr(cbind(global$beta, local$beta[[m]]))$rank, unname(rank$l1 + rank$l2[m]))
  }
})

# Test-only oracle for the original scalar filtering and innovation order.
.scalar_series_reference <- function(n_series, TT, spec, phi, innovation, df, burn, filter_length) {
  if (n_series == 0L) {
    return(matrix(0, TT, 0L))
  }
  n_all <- as.numeric(TT) + burn + 2 * filter_length
  innovations <- matrix(.simulation_innovations(n_all * n_series, innovation, df), n_all, n_series)
  keep <- burn + 2 * filter_length + seq_len(TT)
  if (all(spec$d == 0)) {
    series <- innovations
    if (n_all > 1L) {
      for (t in seq.int(2L, n_all)) {
        series[t, ] <- phi * series[t - 1L, ] + innovations[t, ]
      }
    }
    return(series[keep, , drop = FALSE])
  }
  short <- as.matrix(stats::filter(innovations, phi^(0:filter_length), method = "convolution", sides = 1L))
  result <- matrix(0, TT, n_series)
  end <- cumsum(spec$seg_lens)
  start <- c(1L, end[-length(end)] + 1L)
  for (segment in seq_along(spec$d)) {
    coefficients <- numeric(filter_length + 1)
    coefficients[1L] <- 1
    for (w in seq_len(filter_length)) {
      coefficients[w + 1L] <- coefficients[w] * (w - 1 + spec$d[segment]) / w
    }
    first <- keep[start[segment]]
    last <- keep[end[segment]]
    block <- short[seq.int(first - filter_length, last), , drop = FALSE]
    if (spec$d[segment] != 0 && filter_length > 0L) {
      block <- as.matrix(stats::filter(block, coefficients, method = "convolution", sides = 1L))
    }
    result[seq.int(start[segment], end[segment]), ] <-
      block[seq.int(nrow(block) - spec$seg_lens[segment] + 1L, nrow(block)), , drop = FALSE]
  }
  result
}

test_that("scalar AR(1) simulation agrees exactly with the established engine", {
  args <- .generator_args()
  args$store_components <- TRUE
  for (family in c("Gaussian", "Student-t")) {
    for (nonzero in c(FALSE, TRUE)) {
      args$innovation <- family
      args$memory <- if (nonzero) {
        list(mu_global = c(0.1, 0.2), mu_local = c(-0.2, 0.1), alpha = c(0.2, 0),
             beta = c(0, 0.1), G = c(0.1, 0.2), F = c(0.2, 0.1))
      } else NULL
      actual <- do.call(gen_MMEFM, args)
      expected <- local({
        testthat::local_mocked_bindings(.simulation_series = .scalar_series_reference, .package = "MMEFM")
        do.call(gen_MMEFM, args)
      })
      expect_identical(actual, expected)
    }
  }
})

test_that("stable AR(p) drives latent and error components without changing contracts", {
  args <- .generator_args()
  args$store_components <- TRUE
  args$phi <- c(1.2, -0.5)
  args$error_phi <- c(0.4, -0.2)
  args$memory <- list(mu_global = c(0, 0.1), mu_local = c(-0.2, 0.1), alpha = c(0.2, 0),
                      beta = c(0, 0.1), G = c(0.1, 0.2), F = c(0.2, 0.1))
  set.seed(71L)
  saved <- .Random.seed
  x <- do.call(gen_MMEFM, args)
  expect_identical(.Random.seed, saved)
  expect_identical(x, do.call(gen_MMEFM, args))
  expect_identical(names(x), c("Xt", "dimensions", "rank", "main_effect", "common_component",
                              "error", "memory", "innovation", "components"))
  for (m in seq_len(3L)) {
    expect_true(all(is.finite(x$Xt[[m]])))
    expect_equal(Reduce(`+`, lapply(x$components, function(values) values[[m]])), x$Xt[[m]], tolerance = 1e-13)
  }
  latent_only <- args
  latent_only$error_phi <- 0.2
  expect_identical(do.call(gen_MMEFM, latent_only)$main_effect, x$main_effect)
  expect_identical(do.call(gen_MMEFM, latent_only)$common_component, x$common_component)
  error_only <- args
  error_only$phi <- 0.6
  expect_identical(do.call(gen_MMEFM, error_only)$error, x$error)
  for (name in c("phi", "error_phi")) {
    for (value in list(c(1, 0), c(1.4, 0.1), numeric(), matrix(0.5), c(0.2, NA))) {
      invalid <- args
      invalid[[name]] <- value
      expect_error(do.call(gen_MMEFM, invalid), name)
    }
  }
  args$memory <- NULL
  args$error_rank <- list(ker = 0L, kec = 0L, ker_m = rep(0L, 3L), kec_m = rep(0L, 3L))
  zero_nuisance <- do.call(gen_MMEFM, args)
  expect_identical(dim(zero_nuisance$error$global$G), c(8L, 0L, 0L))
  expect_true(all(is.finite(unlist(zero_nuisance$Xt))))
})

test_that("strength defaults are independent and positional", {
  args <- .generator_args()
  baseline <- do.call(gen_MMEFM, args)
  args$strength <- list()
  expect_identical(do.call(gen_MMEFM, args), baseline)
  args$strength <- list(main_global = NULL, main_local = NULL, common_global = NULL, common_local = NULL)
  expect_identical(do.call(gen_MMEFM, args), baseline)
  args$strength <- list(main_global = c(0.8, 0.9))
  changed <- do.call(gen_MMEFM, args)
  expect_identical(changed$main_effect$local, baseline$main_effect$local)
  expect_identical(changed$common_component, baseline$common_component)
  expect_identical(changed$error, baseline$error)
  for (m in seq_along(args$p)) {
    expect_equal(changed$main_effect$global$A1[[m]], baseline$main_effect$global$A1[[m]] * args$p[m]^(-0.1))
    expect_equal(changed$main_effect$global$B1[[m]], baseline$main_effect$global$B1[[m]] * baseline$dimensions$q[m]^(-0.05))
  }
  args$strength <- list(main_local = list(c(0.8, 1), c(1, 0.9), c(0.9, 0.95)))
  local <- do.call(gen_MMEFM, args)
  expect_identical(local$main_effect$global, baseline$main_effect$global)
  for (m in seq_along(args$p)) {
    xi <- args$strength$main_local[[m]]
    expect_equal(local$main_effect$local$A2[[m]], baseline$main_effect$local$A2[[m]] * args$p[m]^(-(1 - xi[1L]) / 2))
    expect_equal(local$main_effect$local$B2[[m]], baseline$main_effect$local$B2[[m]] * baseline$dimensions$q[m]^(-(1 - xi[2L]) / 2))
  }
  args$q <- args$q[3:1]
  expect_identical(do.call(gen_MMEFM, args), local)
})

test_that("common strengths use group-first row-first column-specific scaling", {
  args <- .generator_args()
  args$p <- c(a = 5L, b = 6L, c = 7L)
  args$q <- c(a = 6L, b = 7L, c = 5L)
  args$rank$kr <- args$rank$kc <- 2L
  args$rank$kr_m <- c(1L, 2L, 1L)
  args$rank$kc_m <- c(2L, 1L, 2L)
  baseline <- do.call(gen_MMEFM, args)
  args$strength <- list(
    common_global = list(list(c(0.7, 0.9), c(1, 0.8)), list(c(0.8, 1), c(0.9, 1)), list(c(0.9, 0.9), c(0.75, 1))),
    common_local = list(list(0.7, c(0.9, 1)), list(c(0.7, 0.95), 0.9), list(1, c(0.8, 0.95))))
  x <- do.call(gen_MMEFM, args)
  expect_identical(x$main_effect, baseline$main_effect)
  expect_identical(x$error, baseline$error)
  for (scope in c("global", "local")) {
    sides <- if (scope == "global") c("Q", "J") else c("R", "C")
    for (m in seq_along(args$p)) {
      for (side in seq_len(2L)) {
        xi <- args$strength[[paste0("common_", scope)]][[m]][[side]]
        loading <- baseline$common_component[[scope]][[sides[side]]][[m]]
        if (length(xi) == 1L) xi <- rep(xi, ncol(loading))
        n <- if (side == 1L) args$p[m] else args$q[m]
        expect_equal(x$common_component[[scope]][[sides[side]]][[m]],
                     sweep(loading, 2L, n^(-(1 - xi) / 2), "*"))
      }
    }
  }
  expanded <- args
  expanded$strength$common_local[[1L]][[1L]] <- rep(0.7, args$rank$kr_m[1L])
  expanded$strength$common_local[[2L]][[2L]] <- rep(0.9, args$rank$kc_m[2L])
  expanded$strength$common_global[[3L]][[1L]] <- 0.9
  expect_identical(do.call(gen_MMEFM, expanded), x)
})

test_that("strength validation rejects names, malformed nesting, and invalid values", {
  args <- .generator_args()
  invalid <- list(1, list(1), list(unknown = 1), list(main_global = c(1, 1), main_global = c(1, 1)),
    setNames(list(c(1, 1)), NA_character_), setNames(list(c(1, 1)), ""),
    list(main_global = c(row = 1, column = 1)), list(main_global = matrix(1, 1L, 2L)),
    list(main_global = 1), list(main_local = list(c(1, 1))),
    list(main_local = structure(rep(list(c(1, 1)), 3L), dim = 3L)),
    list(common_local = rep(list(structure(list(1, 1), dim = 2L)), 3L)),
    list(main_local = list(a = c(1, 1), b = c(1, 1), c = c(1, 1))),
    list(main_local = list(c(row = 1, column = 1), c(1, 1), c(1, 1))),
    list(common_global = list(list(1), list(1, 1), list(1, 1))),
    list(common_local = list(list(1, 1, 1), list(1, 1), list(1, 1))),
    list(common_global = list(list(row = 1, column = 1), list(1, 1), list(1, 1))),
    list(common_global = list(a = list(1, 1), b = list(1, 1), c = list(1, 1))),
    list(common_local = list(list(c(a = 1), 1), list(1, 1), list(1, 1))),
    list(common_global = list(list(c(0.8, 0.9), 1), list(1, 1), list(1, 1))))
  for (strength in invalid) {
    args$strength <- strength
    expect_error(do.call(gen_MMEFM, args), "strength")
  }
  for (value in c(0.5, 0.4, 1.1, NA_real_, Inf, -Inf)) {
    for (field in c("main_global", "main_local", "common_global", "common_local")) {
      entry <- switch(field, main_global = c(value, 1), main_local = rep(list(c(value, 1)), 3L),
                      rep(list(list(value, 1)), 3L))
      args$strength <- setNames(list(entry), field)
      expect_error(do.call(gen_MMEFM, args), paste0("strength\\$", field))
    }
  }
  args$strength <- list(main_global = c(0.5001, 1))
  expect_no_warning(do.call(gen_MMEFM, args))
  args$strength <- list(common_local = list(list(1, numeric()), list(1, 1), list(1, 1)))
  expect_error(do.call(gen_MMEFM, args), "common_local.*1.*column")
})

test_that("error ranks default independently and match canonical group order", {
  args <- .generator_args()
  baseline <- do.call(gen_MMEFM, args)
  args$error_rank <- list()
  expect_identical(do.call(gen_MMEFM, args), baseline)
  args$error_rank <- list(ker = 2L, kec = 2L, ker_m = rep(2L, 3L), kec_m = rep(2L, 3L))
  expect_identical(do.call(gen_MMEFM, args), baseline)
  args$error_rank <- list(ker_m = c(1L, 2L, 0L), kec_m = c(2L, 1L, 0L))
  partial <- do.call(gen_MMEFM, args)
  expect_identical(dim(partial$error$global$G), c(8L, 2L, 2L))
  for (m in seq_along(args$p)) {
    expect_identical(dim(partial$error$local$F[[m]]),
                     c(8L, args$error_rank$ker_m[m], args$error_rank$kec_m[m]))
  }
  args$error_rank <- list(ker = 2, kec = 2, ker_m = c(c = 0, a = 1, b = 2), kec_m = c(b = 1, c = 0, a = 2))
  expect_identical(do.call(gen_MMEFM, args), partial)
  args$error_rank <- list(ker = 1L)
  x <- do.call(gen_MMEFM, args)
  expect_identical(dim(x$error$global$G), c(8L, 1L, 2L))
  expect_identical(dim(x$error$local$F$a), c(8L, 2L, 2L))
  args$error_rank <- list(ker = min(args$p), kec = min(args$q), ker_m = unname(args$p), kec_m = unname(args$q[names(args$p)]))
  expect_no_warning(do.call(gen_MMEFM, args))
  args$p <- unname(args$p)
  args$q <- unname(args$q)
  args$error_rank <- list(ker_m = c(c = 1L, b = 2L, a = 0L))
  expect_identical(dim(do.call(gen_MMEFM, args)$error$local$F[[1L]]), c(8L, 1L, 2L))
})

test_that("error ranks reject malformed or infeasible fields", {
  args <- .generator_args()
  invalid <- list(1, list(1), list(unknown = 1), list(ker = 1, ker = 2),
    setNames(list(1), NA_character_), setNames(list(1), ""), list(ker = NULL),
    list(ker = c(1, 2)), list(ker_m = 1), list(kec_m = rep(1, 4L)),
    list(ker_m = c(a = 1, a = 2, b = 0)), list(kec_m = c(a = 1, b = 1, d = 1)),
    list(ker_m = setNames(c(1, 1, 1), c("a", "b", ""))))
  for (error_rank in invalid) {
    args$error_rank <- error_rank
    expect_error(do.call(gen_MMEFM, args), "error_rank")
  }
  for (field in c("ker", "kec", "ker_m", "kec_m")) {
    for (value in c(-1, 1.5, NA_real_, Inf)) {
      args$error_rank <- setNames(list(if (grepl("_m$", field)) rep(value, 3L) else value), field)
      expect_error(do.call(gen_MMEFM, args), paste0("error_rank\\$", field))
    }
    args$error_rank <- setNames(list(if (grepl("_m$", field)) c(99L, 0L, 0L) else 99L), field)
    expect_error(do.call(gen_MMEFM, args), paste0("error_rank\\$", field, ".*group.*a"))
  }
})
