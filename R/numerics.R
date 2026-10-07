.leading_eigenvectors <- function(S, k) {
  eigen((S + t(S)) / 2, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
}

# Unregularized singular Gram systems remain visible failures.
.solve_gram <- function(A, B, lambda = 0) {
  if (lambda > 0) A <- A + lambda * diag(nrow(A))
  solve(A, B)
}
