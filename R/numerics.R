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
