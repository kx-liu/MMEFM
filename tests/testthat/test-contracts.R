test_that("Xt metadata preserves unequal dimensions and optional group names", {
  Xt <- list(a = array(0, c(1, 3, 4)), b = array(1L, c(1, 5, 2)))
  original <- Xt
  expect_identical(.validate_Xt(Xt), list(M = 2L, T = 1L, p = c(a = 3L, b = 5L),
                                        q = c(a = 4L, b = 2L), group_names = c("a", "b")))
  expect_identical(Xt, original)
  expect_identical(.validate_Xt(unname(Xt)), list(M = 2L, T = 1L, p = c(3L, 5L),
                                                q = c(4L, 2L), group_names = NULL))
})

test_that("Xt rejects noncanonical structures and invalid dimensions", {
  Xt <- list(array(0, c(2, 3, 4)), array(0, c(2, 5, 3)))
  expect_error(.validate_Xt(1:2), "list")
  expect_error(.validate_Xt(Xt[1]), "at least two")
  expect_error(.validate_Xt(list()), "at least two")
  for (bad in list(matrix(0, 2, 3), data.frame(x = 1:2), 1:2,
                   array("x", c(2, 3, 4)), array(TRUE, c(2, 3, 4)))) {
    invalid <- Xt
    invalid[[1]] <- bad
    expect_error(.validate_Xt(invalid), "numeric three-dimensional")
  }
  invalid <- Xt
  invalid[[2]] <- array(0, c(3, 5, 3))
  expect_error(.validate_Xt(invalid), "same time dimension")
  for (dims in list(c(0, 3, 4), c(2, 1, 4), c(2, 3, 1))) {
    invalid <- Xt
    invalid[[1]] <- array(0, dims)
    expect_error(.validate_Xt(invalid), "T >= 1, p_m >= 2, and q_m >= 2")
  }
})

test_that("Xt requires finite entries and unambiguous group names", {
  Xt <- list(a = array(0, c(1, 2, 2)), b = array(0, c(1, 2, 2)))
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    invalid <- Xt
    invalid[[1]][1] <- value
    expect_error(.validate_Xt(invalid), "finite entries")
  }
  for (group_names in list(c("a", "a"), c("a", ""), c("a", NA_character_))) {
    invalid <- Xt
    names(invalid) <- group_names
    expect_error(.validate_Xt(invalid), "nonempty, nonmissing, and unique")
  }
})

test_that("complete ranks are canonical integers in Xt group order", {
  data_info <- .validate_Xt(list(a = array(0, c(1, 3, 5)), b = array(0, c(1, 5, 4))))
  rank <- list(r1 = 1, l1 = 3, r2 = c(2, 4), l2 = c(4, 3),
               kr = 2, kc = 1, kr_m = c(1, 3), kc_m = c(2, 1))
  expected <- list(r1 = 1L, l1 = 3L, r2 = c(a = 2L, b = 4L), l2 = c(a = 4L, b = 3L),
                   kr = 2L, kc = 1L, kr_m = c(a = 1L, b = 3L), kc_m = c(a = 2L, b = 1L))
  expect_identical(.validate_rank(rank, data_info), expected)
  for (component in c("r2", "l2", "kr_m", "kc_m")) names(rank[[component]]) <- c("a", "b")
  expect_identical(.validate_rank(rank, data_info), expected)
  for (component in c("r2", "l2", "kr_m", "kc_m")) rank[[component]] <- rank[[component]][c("b", "a")]
  expect_identical(.validate_rank(rank, data_info), expected)
  expect_identical(.validate_rank(rank[rev(names(rank))], data_info), expected)

  unnamed_info <- .validate_Xt(list(array(0, c(1, 3, 5)), array(0, c(1, 5, 4))))
  # Names cannot override position when Xt itself has no names.
  positional <- list(r1 = 1L, l1 = 1L, r2 = c(b = 2L, a = 4L), l2 = c(b = 3L, a = 2L),
                     kr = 1L, kc = 1L, kr_m = c(b = 2L, a = 4L), kc_m = c(b = 3L, a = 2L))
  canonical <- .validate_rank(positional, unnamed_info)
  expect_identical(canonical$r2, c(2L, 4L))
  expect_identical(canonical$l2, c(3L, 2L))
  expect_null(names(canonical$kr_m))
  expect_null(names(canonical$kc_m))
})

test_that("manual ranks require exactly the eight canonical components", {
  data_info <- .validate_Xt(list(array(0, c(1, 3, 5)), array(0, c(1, 5, 4))))
  rank <- list(r1 = 1L, l1 = 1L, r2 = c(1L, 1L), l2 = c(1L, 1L),
               kr = 1L, kc = 1L, kr_m = c(1L, 1L), kc_m = c(1L, 1L))
  expect_identical(.validate_rank(rank, data_info), rank)
  expect_error(.validate_rank(NULL, data_info), "complete named list")
  expect_error(.validate_rank(unname(rank), data_info), "complete named list")
  expect_error(.validate_rank(rank[-1], data_info), "complete named list")
  expect_error(.validate_rank(rank[1:2], data_info), "complete named list")
  expect_error(.validate_rank(c(rank, list(other = 1)), data_info), "complete named list")
  for (alias in c("k_global", "k_local", "rl_global", "rl_local", "k_r", "k_c", "k_r_m", "k_c_m")) {
    invalid <- rank
    names(invalid)[5] <- alias
    expect_error(.validate_rank(invalid, data_info), "complete named list")
  }
  for (bad_names in list(rep("r1", 8), c(names(rank)[-8], NA_character_))) {
    invalid <- rank
    names(invalid) <- bad_names
    expect_error(.validate_rank(invalid, data_info), "complete named list")
  }
})

test_that("every rank must be a positive whole-number vector of the right length", {
  data_info <- .validate_Xt(list(array(0, c(1, 3, 5)), array(0, c(1, 5, 4))))
  rank <- list(r1 = 1L, l1 = 1L, r2 = c(1L, 1L), l2 = c(1L, 1L),
               kr = 1L, kc = 1L, kr_m = c(1L, 1L), kc_m = c(1L, 1L))
  for (component in names(rank)) {
    for (value in c(1.5, 0, -1, NA_real_, NaN, Inf, -Inf)) {
      invalid <- rank
      invalid[[component]][1] <- value
      expect_error(.validate_rank(invalid, data_info), "positive whole-number")
    }
    invalid <- rank
    invalid[[component]] <- rep(1, length(rank[[component]]) + 1L)
    expect_error(.validate_rank(invalid, data_info), "positive whole-number")
    invalid[[component]] <- matrix(rank[[component]], nrow = 1)
    expect_error(.validate_rank(invalid, data_info), "positive whole-number")
    invalid[[component]] <- as.character(rank[[component]])
    expect_error(.validate_rank(invalid, data_info), "positive whole-number")
  }
})

test_that("global and local rank bounds follow their respective dimensions", {
  data_info <- .validate_Xt(list(array(0, c(1, 3, 5)), array(0, c(1, 5, 4))))
  rank <- list(r1 = 1L, l1 = 1L, r2 = c(1L, 1L), l2 = c(1L, 1L),
               kr = 1L, kc = 1L, kr_m = c(1L, 1L), kc_m = c(1L, 1L))
  for (component in c("r1", "l1", "kr", "kc")) {
    limit <- if (component %in% c("r1", "kr")) min(data_info$p) else min(data_info$q)
    for (value in c(limit, limit + 1L)) {
      invalid <- rank
      invalid[[component]] <- value
      expect_error(.validate_rank(invalid, data_info), "minimum group dimension")
    }
  }
  for (component in c("r2", "l2", "kr_m", "kc_m")) {
    dimensions <- if (component %in% c("r2", "kr_m")) data_info$p else data_info$q
    for (m in seq_len(data_info$M)) {
      for (value in c(dimensions[m], dimensions[m] + 1L)) {
        invalid <- rank
        invalid[[component]][m] <- value
        expect_error(.validate_rank(invalid, data_info), "corresponding group dimension")
      }
    }
  }
})

test_that("named group ranks must be complete and unambiguous", {
  data_info <- .validate_Xt(list(a = array(0, c(1, 3, 5)), b = array(0, c(1, 5, 4))))
  rank <- list(r1 = 1L, l1 = 1L, r2 = c(1L, 1L), l2 = c(1L, 1L),
               kr = 1L, kc = 1L, kr_m = c(1L, 1L), kc_m = c(1L, 1L))
  for (component in c("r2", "l2", "kr_m", "kc_m")) {
    invalid <- rank
    names(invalid[[component]]) <- c("a", "extra")
    expect_error(.validate_rank(invalid, data_info), "match the Xt group names exactly")
    for (value_names in list(c("a", "a"), c("a", ""), c("a", NA_character_))) {
      names(invalid[[component]]) <- value_names
      expect_error(.validate_rank(invalid, data_info), "nonempty, nonmissing, and unique")
    }
    invalid[[component]] <- c(a = 1L)
    expect_error(.validate_rank(invalid, data_info), "positive whole-number")
  }
})

test_that("the internal rank constructor preserves ranks and optional diagnostics", {
  ranks <- list(r1 = 1L, l1 = 1L, r2 = c(a = 1L, b = 2L), l2 = c(a = 2L, b = 1L),
                kr = 1L, kc = 1L, kr_m = c(a = 1L, b = 2L), kc_m = c(a = 2L, b = 1L))
  object <- .new_mmefm_rank(ranks)
  expect_identical(names(object), c(names(ranks), "diagnostics"))
  expect_identical(class(object), "mmefm_rank")
  expect_identical(unclass(object)[names(ranks)], ranks)
  expect_null(object$diagnostics)
  diagnostics <- list(ratios = c(0.1, 0.9))
  supplied <- .new_mmefm_rank(ranks, diagnostics)
  expect_identical(supplied$diagnostics, diagnostics)
  expect_identical(unclass(supplied)[names(ranks)], ranks)
})
