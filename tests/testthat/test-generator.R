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
  args$error_rank_global <- c(0L, 2L)
  args$error_rank_local <- list(c(1L, 0L), c(0L, 0L), c(0L, 2L))
  zero <- do.call(gen_MMEFM, args)
  expect_identical(dim(zero$error$global$G), c(8L, 0L, 2L))
  expect_identical(dim(zero$error$global$Q$a), c(4L, 0L))
  for (m in seq_len(3L)) {
    expect_true(all(zero$components$error_global[[m]] == 0))
    expect_true(all(zero$components$error_local[[m]] == 0))
    expect_identical(dim(zero$error$local$F[[m]]), c(8L, args$error_rank_local[[m]]))
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
  args$common_strength_global <- list(row = list(c = 0.8, b = c(0.7, 0.9), a = 0.8), column = rep(list(1), 3L))
  scalar <- do.call(gen_MMEFM, args)
  args$common_strength_global$row <- list(c = c(0.8, 0.8), b = c(0.7, 0.9), a = c(0.8, 0.8))
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
  invalid <- list(T = 1.5, p = c(1, 4, 5), q = c(4, 5), main_strength_global = c(0.5, 1),
                  main_strength_local = list(c(1, 2)), common_strength_global = list(row = rep(list(0.4), 3L), column = rep(list(1), 3L)),
                  common_strength_local = list(row = list(1)), phi = 1, error_phi = Inf,
                  error_rank_global = c(-1, 2), error_rank_local = rep(list(c(1.5, 2)), 3L),
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
