.leading_eigenvectors <- function(S, k) {
  eigen((S + t(S)) / 2, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
}

# Unregularized singular Gram systems remain visible failures.
.solve_gram <- function(A, B, lambda = 0) {
  if (lambda > 0) A <- A + lambda * diag(nrow(A))
  solve(A, B)
}

# A and B are equal-rank orthonormal loading bases.
.subspace_distance <- function(A, B) {
  singular_values <- svd(crossprod(A, B), nu = 0L, nv = 0L)$d
  min_cos <- min(pmax(0, pmin(1, singular_values)))
  sqrt(1 - min_cos^2)
}

# Shared stochastic-boundary policy; low-level algorithms consume the caller's stream.
.with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}
