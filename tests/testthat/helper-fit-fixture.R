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
