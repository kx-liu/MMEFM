# Group order is shared by strengths, memory specifications, and nuisance ranks.
.simulation_groups <- function(x, M, group_names, label) {
  if (!is.list(x) || is.data.frame(x) || length(x) != M) stop(label, " must be a list of length M.")
  if (!is.null(names(x))) {
    if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) stop(label, " group names must be nonempty and unique.")
    if (!is.null(group_names)) {
      if (!setequal(names(x), group_names)) stop(label, " names must match the dimension group names.")
      x <- x[match(group_names, names(x))]
    }
  }
  names(x) <- group_names
  x
}

.simulation_memory_spec <- function(spec, TT, label) {
  if (is.numeric(spec) && is.null(dim(spec))) {
    d <- spec
    if (!length(d) || length(d) > TT) stop(label, " must have between one and T segments.")
    seg_lens <- rep(as.integer(floor(TT / length(d))), length(d))
    seg_lens[length(d)] <- TT - sum(seg_lens[-length(d)])
  } else if (is.list(spec) && length(spec) == 2L && !is.null(names(spec)) &&
             !anyNA(names(spec)) && !anyDuplicated(names(spec)) && setequal(names(spec), c("d", "seg_lens"))) {
    d <- spec$d
    seg_lens <- spec$seg_lens
  } else stop(label, " must be numeric memory parameters or list(d, seg_lens).")
  if (!is.numeric(d) || !is.null(dim(d)) || !length(d) || any(!is.finite(d)) || any(d >= 0.5)) stop(label, "$d must contain finite memory parameters below 0.5.")
  if (!is.numeric(seg_lens) || !is.null(dim(seg_lens)) || length(seg_lens) != length(d) ||
      any(!is.finite(seg_lens)) || any(seg_lens < 1) || any(seg_lens != floor(seg_lens)) || sum(seg_lens) != TT) {
    stop(label, "$seg_lens must be positive whole numbers matching d and summing to T.")
  }
  list(d = as.numeric(d), seg_lens = as.integer(seg_lens))
}

.simulation_memory <- function(memory, TT, M, group_names) {
  fields <- c("mu_global", "mu_local", "alpha", "beta", "G", "F")
  if (is.null(memory)) memory <- stats::setNames(rep(list(c(0, 0)), length(fields)), fields)
  if (!is.list(memory) || length(memory) != length(fields) || is.null(names(memory)) ||
      anyNA(names(memory)) || anyDuplicated(names(memory)) || !setequal(names(memory), fields)) {
    stop("memory must contain exactly mu_global, mu_local, alpha, beta, G, and F.")
  }
  memory <- memory[fields]
  for (name in c("mu_global", "mu_local", "G")) memory[[name]] <- .simulation_memory_spec(memory[[name]], TT, paste0("memory$", name))
  for (name in c("alpha", "beta", "F")) {
    spec <- memory[[name]]
    split <- name != "F" && is.list(spec) && !is.null(names(spec)) && any(names(spec) %in% c("global", "local"))
    if (split) {
      if (length(spec) != 2L || anyNA(names(spec)) || anyDuplicated(names(spec)) || !setequal(names(spec), c("global", "local"))) stop("memory$", name, " must contain exactly global and local.")
      global <- .simulation_memory_spec(spec$global, TT, paste0("memory$", name, "$global"))
      local <- spec$local
    } else {
      local <- spec
      if (name != "F") global <- .simulation_memory_spec(spec, TT, paste0("memory$", name))
    }
    basic <- is.numeric(local) || (is.list(local) && !is.null(names(local)) && setequal(names(local), c("d", "seg_lens")))
    if (basic) {
      local <- rep(list(.simulation_memory_spec(local, TT, paste0("memory$", name, "$local"))), M)
      names(local) <- group_names
    } else {
      local <- .simulation_groups(local, M, group_names, paste0("memory$", name, "$local"))
      for (m in seq_len(M)) local[[m]] <- .simulation_memory_spec(local[[m]], TT, paste0("memory$", name, "$local[[", m, "]]"))
    }
    memory[[name]] <- if (name == "F") local else list(global = global, local = local)
  }
  memory
}

.simulation_innovations <- function(n, innovation, df) {
  if (innovation == "Gaussian") stats::rnorm(n) else stats::rt(n, df = df)
}

# The retained sample follows two full filter warmups, as in the numerical DGP.
.simulation_series <- function(n_series, TT, spec, phi, innovation, df, burn, filter_length) {
  if (n_series == 0L) return(matrix(0, TT, 0L))
  n_all <- as.numeric(TT) + burn + 2 * filter_length
  innovations <- matrix(.simulation_innovations(n_all * n_series, innovation, df), n_all, n_series)
  keep <- burn + 2 * filter_length + seq_len(TT)
  if (all(spec$d == 0)) {
    series <- innovations
    if (n_all > 1L) for (t in seq.int(2L, n_all)) series[t, ] <- phi * series[t - 1L, ] + innovations[t, ]
    return(series[keep, , drop = FALSE])
  }
  short <- as.matrix(stats::filter(innovations, phi^(0:filter_length), method = "convolution", sides = 1L))
  result <- matrix(0, TT, n_series)
  end <- cumsum(spec$seg_lens)
  start <- c(1L, end[-length(end)] + 1L)
  for (segment in seq_along(spec$d)) {
    coefficients <- numeric(filter_length + 1)
    coefficients[1L] <- 1
    if (filter_length > 0L) {
      for (w in seq_len(filter_length)) coefficients[w + 1L] <- coefficients[w] * (w - 1 + spec$d[segment]) / w
    }
    first <- keep[start[segment]]
    last <- keep[end[segment]]
    block <- short[seq.int(first - filter_length, last), , drop = FALSE]
    if (spec$d[segment] != 0 && filter_length > 0L) block <- as.matrix(stats::filter(block, coefficients, method = "convolution", sides = 1L))
    result[seq.int(start[segment], end[segment]), ] <- block[seq.int(nrow(block) - spec$seg_lens[segment] + 1L, nrow(block)), , drop = FALSE]
  }
  result
}

.simulation_loading <- function(n, k, strength) {
  if (length(strength) == 1L) strength <- rep(strength, k)
  loading <- matrix(stats::rnorm(as.numeric(n) * k), n, k)
  loading <- sweep(loading, 2L, n^(-(1 - strength) / 2), "*")
  sweep(loading, 2L, colMeans(loading), "-")
}

.simulation_error_loading <- function(n, k, zero_prob) {
  loading <- matrix(stats::rnorm(as.numeric(n) * k), n, k)
  mask <- matrix(stats::runif(as.numeric(n) * k) < zero_prob, n, k)
  loading[mask] <- 0
  loading
}

#' Simulate a Multilevel Main Effects Matrix Factor Model
#'
#' Generate canonical group arrays and nested statistical truth using the manuscript's finite-sample numerical DGP.
#' @param T Positive whole-number time dimension within R's integer range.
#' @param p,q Whole-number row and column dimension vectors, each at least two, with the same length M >= 2. Both must be unnamed or have unique nonempty matching group names; q is reordered to p order.
#' @param rank Complete canonical eight-field rank list or `mmefm_rank`; see [est_MMEFM()]. Positive nominal model ranks apply. All model loadings are centered under IC1, so `r1 + r2[m]` and `kr + kr_m[m]` must be at most `p_m - 1`, and `l1 + l2[m]` and `kc + kc_m[m]` at most `q_m - 1`; equality at these centered dimensions is allowed. Supplied diagnostics are discarded. Each global/local main-score rank sum must also be at most T; equality is allowed.
#' @param main_strength_global Length-two row/column main-loading strengths in (0.5, 1].
#' @param main_strength_local NULL to repeat global strengths, or a group list of length-two strengths.
#' @param common_strength_global,common_strength_local NULL for unit strengths, or exactly `list(row = ..., column = ...)`, each a group list. Entries are scalar or one strength per corresponding loading column, in (0.5, 1]. Named group lists are matched to p order; unnamed lists are positional.
#' @param memory NULL for two zero-memory segments for all latent components, or exactly a list with `mu_global`, `mu_local`, `alpha`, `beta`, `G`, `F`. A basic spec is a numeric vector of finite d < 0.5 or `list(d = ..., seg_lens = ...)` with positive whole segment lengths summing to T. Numeric specs split T equally with remainder in the final segment. Means and G use basic specs. F uses a basic spec or group list of basic specs. Alpha/beta use a basic spec or `list(global = <basic>, local = <basic or group list>)`. Local mean innovations share one spec because group demeaning mixes them.
#' @param phi,error_phi Finite geometric short-memory coefficients with absolute value below one, for latent and error series respectively.
#' @param error_rank_global Length-two nonnegative whole row/column nuisance ranks fitting all groups.
#' @param error_rank_local NULL for c(2, 2) in every group, or a group list of nonnegative whole row/column nuisance ranks fitting each group. Zero nuisance ranks produce zero components.
#' @param error_loading_zero_prob Probability in [0, 1) of independently zeroing each Gaussian error-loading entry.
#' @param innovation Gaussian (default) or raw Student-t temporal innovations. Loadings and error scale entries always use Gaussian draws.
#' @param df Finite positive Student-t degrees of freedom, retained in metadata even for Gaussian. Raw Student-t draws are not variance-standardized; theoretical fourth-moment conditions require an appropriate df, but simulation permits any positive df.
#' @param burn Nonnegative whole-number burn-in within integer range.
#' @param filter_length NULL for `min(2000, max(200, 5*T))`, or a nonnegative whole fractional/geometric filter truncation length. Zero-memory series use AR(1) recursion; other series use geometric convolution followed by segment-specific fractional filters.
#' @param active_global NULL for all groups, a nonmissing logical length-M vector, or unique whole group indices (empty indices allowed). Zero or at least two active groups are permitted; a singleton errors. Ownership is applied after ordinary global draws by zeroing inactive Q/J. Nominal positive ranks and G remain unchanged.
#' @param seed Nonnegative whole-number seed within integer range. Caller RNG state is preserved on success and error; RNG kind is unchanged.
#' @param store_components One nonmissing logical. FALSE avoids storing seven additional full-array component lists; truth loadings and factor series are always returned.
#' @details Loading column entries are independent Gaussian draws scaled by n^(-(1-xi)/2), then demeaned for IC1. Spatial global/local loadings are not artificially orthogonalized. Local grand means sum to zero across groups at each time. Local alpha/beta score matrices are projected with ordinary solve onto the complement of the global temporal score span, enforcing exact sample IC2 orthogonality to roundoff; singular Grams error without ridge or fallback.
#'
#' Error loadings are independent sparse Gaussian matrices across groups. The numerical DGP deliberately relaxes E1 cross-group error-loading orthogonality, corresponding to the manuscript's relaxed specification. Error AR processes are not rescaled to unit variance, and Sigma_epsilon entries are raw absolute Gaussian draws. This finite-sample simulator does not literally enforce every asymptotic E1/E2 regularity condition.
#' @return An ordinary list with exactly `Xt`, `dimensions`, `rank`, `main_effect`, `common_component`, `error`, `memory`, `innovation`, and `components`. Dimensions contain M, T, p, q and optional group names. Rank is an `mmefm_rank` with NULL diagnostics. Main truth contains `direct` c/a/b, `global` mu/A1/alpha/B1/beta, and `local` mu/A2/alpha/B2/beta. Common truth contains global Q/J/G/active and local R/C/F. Error truth contains global Q/J/G, local R/C/F, and group `sigma` matrices. Scores and means are time-by-direction matrices; cores are time-by-row-rank-by-column-rank arrays. All group lists preserve optional names.
#'
#' Memory contains canonical d/seg_lens specs, with alpha/beta split into global and group-local specs and F expanded by group. Innovation contains exactly family and df. Components is NULL or exactly seven group-array lists: global_main, local_main, global_common, local_common, error_global, error_local, error_idiosyncratic; their sum recovers Xt.
#' @examples
#' rank <- list(r1 = 1, l1 = 1, r2 = c(1, 1), l2 = c(1, 1),
#'              kr = 1, kc = 1, kr_m = c(1, 1), kc_m = c(1, 1))
#' simulation <- gen_MMEFM(8, c(4, 5), c(5, 4), rank,
#'                         burn = 4, filter_length = 8, store_components = FALSE)
#' dim(simulation$Xt[[1]])
#' @export
gen_MMEFM <- function(T, p, q, rank, main_strength_global = c(1, 1), main_strength_local = NULL,
                      common_strength_global = NULL, common_strength_local = NULL, memory = NULL,
                      phi = 0.6, error_phi = 0.2, error_rank_global = c(2L, 2L), error_rank_local = NULL,
                      error_loading_zero_prob = 0.95, innovation = c("Gaussian", "Student-t"), df = 6,
                      burn = 200L, filter_length = NULL, active_global = NULL, seed = 2026L, store_components = FALSE) {
  innovation <- match.arg(innovation)
  controls <- list(T = T, burn = burn, seed = seed)
  if (!is.null(filter_length)) controls$filter_length <- filter_length
  for (name in names(controls)) {
    value <- controls[[name]]
    minimum <- if (name == "T") 1L else 0L
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != floor(value) || value < minimum || value > .Machine$integer.max) stop(name, " must be one finite whole number within integer range, >= ", minimum, ".")
  }
  T <- as.integer(T)
  burn <- as.integer(burn)
  seed <- as.integer(seed)
  if (is.null(filter_length)) filter_length <- min(2000, max(200, 5 * as.numeric(T)))
  filter_length <- as.integer(filter_length)
  for (name in c("p", "q")) {
    value <- get(name)
    if (!is.numeric(value) || !is.null(dim(value)) || length(value) < 2L || any(!is.finite(value)) || any(value < 2) || any(value != floor(value)) || any(value > .Machine$integer.max)) stop(name, " must be a vector of at least two whole dimensions >= 2 within integer range.")
  }
  if (length(p) != length(q)) stop("p and q must have the same length M.")
  group_names <- names(p)
  if (is.null(group_names) != is.null(names(q))) stop("p and q must both be named or both unnamed.")
  if (!is.null(group_names)) {
    for (name in list(group_names, names(q))) if (anyNA(name) || any(!nzchar(name)) || anyDuplicated(name)) stop("p/q group names must be nonempty, nonmissing, and unique.")
    if (!setequal(group_names, names(q))) stop("p/q group names must match exactly.")
    q <- q[match(group_names, names(q))]
  }
  p <- as.integer(p)
  q <- as.integer(q)
  names(p) <- names(q) <- group_names
  M <- length(p)
  dimensions <- list(M = M, T = T, p = p, q = q, group_names = group_names)
  ranks <- .validate_rank(rank, dimensions)
  for (side in c("row", "column")) {
    total <- if (side == "row") as.numeric(ranks$r1) + ranks$r2 else as.numeric(ranks$l1) + ranks$l2
    if (any(total > T)) stop("Global/local ", side, " main-score rank sum must not exceed T in group(s): ", paste(if (is.null(group_names)) which(total > T) else group_names[total > T], collapse = ", "), ".")
  }
  if (is.null(main_strength_local)) main_strength_local <- rep(list(main_strength_global), M)
  main_strength_local <- .simulation_groups(main_strength_local, M, group_names, "main_strength_local")
  for (value in c(list(main_strength_global), main_strength_local)) {
    if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 2L || any(!is.finite(value)) || any(value <= 0.5 | value > 1)) stop("main_strength_global/local must have two finite strengths in (0.5, 1].")
  }
  strengths <- list(global = common_strength_global, local = common_strength_local)
  for (scope in names(strengths)) {
    strength <- strengths[[scope]]
    if (is.null(strength)) strength <- list(row = rep(list(1), M), column = rep(list(1), M))
    if (!is.list(strength) || length(strength) != 2L || is.null(names(strength)) || anyNA(names(strength)) || anyDuplicated(names(strength)) || !setequal(names(strength), c("row", "column"))) stop("common_strength_", scope, " must contain exactly row and column.")
    strength <- strength[c("row", "column")]
    for (side in names(strength)) {
      strength[[side]] <- .simulation_groups(strength[[side]], M, group_names, paste0("common_strength_", scope, "$", side))
      required <- if (scope == "global") rep(if (side == "row") ranks$kr else ranks$kc, M) else if (side == "row") ranks$kr_m else ranks$kc_m
      for (m in seq_len(M)) {
        value <- strength[[side]][[m]]
        if (!is.numeric(value) || !is.null(dim(value)) || !(length(value) %in% c(1L, required[m])) || any(!is.finite(value)) || any(value <= 0.5 | value > 1)) stop("common_strength_", scope, "$", side, "[[", m, "]] must be scalar or match its rank, with strengths in (0.5, 1].")
        if (length(value) == 1L) value <- rep(value, required[m])
        strength[[side]][[m]] <- value
      }
    }
    strengths[[scope]] <- strength
  }
  memory <- .simulation_memory(memory, T, M, group_names)
  for (name in c("phi", "error_phi")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || abs(value) >= 1) stop(name, " must be one finite number with absolute value below one.")
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) stop("df must be one finite positive number.")
  if (!is.numeric(error_loading_zero_prob) || length(error_loading_zero_prob) != 1L || !is.finite(error_loading_zero_prob) || error_loading_zero_prob < 0 || error_loading_zero_prob >= 1) stop("error_loading_zero_prob must be in [0, 1).")
  if (is.null(error_rank_local)) error_rank_local <- rep(list(c(2L, 2L)), M)
  error_rank_local <- .simulation_groups(error_rank_local, M, group_names, "error_rank_local")
  for (m in seq_len(M)) {
    for (scope in c("global", "local")) {
      value <- if (scope == "global") error_rank_global else error_rank_local[[m]]
      if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 2L || any(!is.finite(value)) || any(value < 0) || any(value != floor(value)) || any(value > c(p[m], q[m]))) stop("error_rank_", scope, " must have two nonnegative whole ranks fitting group ", m, ".")
      if (scope == "local") error_rank_local[[m]] <- as.integer(value)
    }
  }
  error_rank_global <- as.integer(error_rank_global)
  active <- rep(TRUE, M)
  if (!is.null(active_global)) {
    if (is.logical(active_global)) {
      if (!is.null(dim(active_global)) || length(active_global) != M || anyNA(active_global)) stop("active_global must be a nonmissing logical vector of length M.")
      active <- unname(active_global)
    } else {
      if (!is.numeric(active_global) || !is.null(dim(active_global)) || any(!is.finite(active_global)) || any(active_global != floor(active_global)) || any(active_global < 1 | active_global > M) || anyDuplicated(active_global)) stop("active_global must contain unique whole group indices between 1 and M.")
      active <- seq_len(M) %in% sort(as.integer(active_global))
    }
  }
  if (sum(active) == 1L) stop("active_global cannot contain a singleton: global factors require at least two groups.")
  names(active) <- group_names
  if (!is.logical(store_components) || length(store_components) != 1L || is.na(store_components)) stop("store_components must be one nonmissing logical.")
  .with_preserved_seed(seed, {
    mu <- .simulation_series(1L, T, memory$mu_global, phi, innovation, df, burn, filter_length)
    local_mu <- .simulation_series(M, T, memory$mu_local, phi, innovation, df, burn, filter_length)
    local_mu <- sweep(local_mu, 1L, rowMeans(local_mu), "-")
    c <- mu_local <- vector("list", M)
    names(c) <- names(mu_local) <- group_names
    for (m in seq_len(M)) {
      mu_local[[m]] <- local_mu[, m, drop = FALSE]
      c[[m]] <- mu + mu_local[[m]]
    }
    A1 <- B1 <- A2 <- B2 <- vector("list", M)
    names(A1) <- names(B1) <- names(A2) <- names(B2) <- group_names
    for (m in seq_len(M)) A1[[m]] <- .simulation_loading(p[m], ranks$r1, main_strength_global[1L])
    for (m in seq_len(M)) B1[[m]] <- .simulation_loading(q[m], ranks$l1, main_strength_global[2L])
    for (m in seq_len(M)) A2[[m]] <- .simulation_loading(p[m], ranks$r2[m], main_strength_local[[m]][1L])
    for (m in seq_len(M)) B2[[m]] <- .simulation_loading(q[m], ranks$l2[m], main_strength_local[[m]][2L])
    alpha <- .simulation_series(ranks$r1, T, memory$alpha$global, phi, innovation, df, burn, filter_length)
    beta <- .simulation_series(ranks$l1, T, memory$beta$global, phi, innovation, df, burn, filter_length)
    alpha_local <- beta_local <- a <- b <- vector("list", M)
    names(alpha_local) <- names(beta_local) <- names(a) <- names(b) <- group_names
    for (m in seq_len(M)) {
      alpha_local[[m]] <- .simulation_series(ranks$r2[m], T, memory$alpha$local[[m]], phi, innovation, df, burn, filter_length)
      beta_local[[m]] <- .simulation_series(ranks$l2[m], T, memory$beta$local[[m]], phi, innovation, df, burn, filter_length)
      alpha_local[[m]] <- alpha_local[[m]] - alpha %*% solve(crossprod(alpha), crossprod(alpha, alpha_local[[m]]))
      beta_local[[m]] <- beta_local[[m]] - beta %*% solve(crossprod(beta), crossprod(beta, beta_local[[m]]))
      a[[m]] <- alpha %*% t(A1[[m]]) + alpha_local[[m]] %*% t(A2[[m]])
      b[[m]] <- beta %*% t(B1[[m]]) + beta_local[[m]] %*% t(B2[[m]])
    }
    G <- array(.simulation_series(as.numeric(ranks$kr) * ranks$kc, T, memory$G, phi, innovation, df, burn, filter_length), c(T, ranks$kr, ranks$kc))
    F <- vector("list", M)
    names(F) <- group_names
    for (m in seq_len(M)) F[[m]] <- array(.simulation_series(as.numeric(ranks$kr_m[m]) * ranks$kc_m[m], T, memory$F[[m]], phi, innovation, df, burn, filter_length), c(T, unname(ranks$kr_m[m]), unname(ranks$kc_m[m])))
    Q <- J <- R <- C <- vector("list", M)
    names(Q) <- names(J) <- names(R) <- names(C) <- group_names
    for (m in seq_len(M)) Q[[m]] <- .simulation_loading(p[m], ranks$kr, strengths$global$row[[m]])
    for (m in seq_len(M)) J[[m]] <- .simulation_loading(q[m], ranks$kc, strengths$global$column[[m]])
    for (m in seq_len(M)) R[[m]] <- .simulation_loading(p[m], ranks$kr_m[m], strengths$local$row[[m]])
    for (m in seq_len(M)) C[[m]] <- .simulation_loading(q[m], ranks$kc_m[m], strengths$local$column[[m]])
    Q_e <- J_e <- R_e <- C_e <- vector("list", M)
    names(Q_e) <- names(J_e) <- names(R_e) <- names(C_e) <- group_names
    for (m in seq_len(M)) Q_e[[m]] <- .simulation_error_loading(p[m], error_rank_global[1L], error_loading_zero_prob)
    for (m in seq_len(M)) J_e[[m]] <- .simulation_error_loading(q[m], error_rank_global[2L], error_loading_zero_prob)
    for (m in seq_len(M)) R_e[[m]] <- .simulation_error_loading(p[m], error_rank_local[[m]][1L], error_loading_zero_prob)
    for (m in seq_len(M)) C_e[[m]] <- .simulation_error_loading(q[m], error_rank_local[[m]][2L], error_loading_zero_prob)
    short_memory <- list(d = 0, seg_lens = T)
    G_e <- array(.simulation_series(prod(error_rank_global), T, short_memory, error_phi, innovation, df, burn, filter_length), c(T, error_rank_global))
    F_e <- sigma <- vector("list", M)
    names(F_e) <- names(sigma) <- group_names
    for (m in seq_len(M)) {
      F_e[[m]] <- array(.simulation_series(prod(error_rank_local[[m]]), T, short_memory, error_phi, innovation, df, burn, filter_length), c(T, error_rank_local[[m]]))
      sigma[[m]] <- matrix(abs(stats::rnorm(as.numeric(p[m]) * q[m])), p[m], q[m])
    }
    for (m in which(!active)) {
      Q[[m]][] <- 0
      J[[m]][] <- 0
    }
    component_names <- c("global_main", "local_main", "global_common", "local_common", "error_global", "error_local", "error_idiosyncratic")
    components <- if (store_components) stats::setNames(vector("list", length(component_names)), component_names) else NULL
    if (store_components) for (name in component_names) {
      components[[name]] <- vector("list", M)
      names(components[[name]]) <- group_names
    }
    Xt <- vector("list", M)
    names(Xt) <- group_names
    for (m in seq_len(M)) {
      pm <- unname(p[m])
      qm <- unname(q[m])
      epsilon <- array(.simulation_series(as.numeric(pm) * qm, T, short_memory, error_phi, innovation, df, burn, filter_length), c(T, pm, qm))
      Xt[[m]] <- array(0, c(T, pm, qm))
      if (store_components) for (name in component_names) components[[name]][[m]] <- array(0, c(T, pm, qm))
      for (t in seq_len(T)) {
        global_main <- matrix(mu[t, 1L], pm, qm) +
          outer(as.vector(A1[[m]] %*% alpha[t, ]), rep(1, qm)) + outer(rep(1, pm), as.vector(B1[[m]] %*% beta[t, ]))
        local_main <- matrix(mu_local[[m]][t, 1L], pm, qm) +
          outer(as.vector(A2[[m]] %*% alpha_local[[m]][t, ]), rep(1, qm)) + outer(rep(1, pm), as.vector(B2[[m]] %*% beta_local[[m]][t, ]))
        global_common <- Q[[m]] %*% matrix(G[t, , , drop = FALSE], ranks$kr, ranks$kc) %*% t(J[[m]])
        local_common <- R[[m]] %*% matrix(F[[m]][t, , , drop = FALSE], ranks$kr_m[m], ranks$kc_m[m]) %*% t(C[[m]])
        error_global <- Q_e[[m]] %*% matrix(G_e[t, , , drop = FALSE], error_rank_global[1L], error_rank_global[2L]) %*% t(J_e[[m]])
        error_local <- R_e[[m]] %*% matrix(F_e[[m]][t, , , drop = FALSE], error_rank_local[[m]][1L], error_rank_local[[m]][2L]) %*% t(C_e[[m]])
        error_idiosyncratic <- sigma[[m]] * matrix(epsilon[t, , , drop = FALSE], pm, qm)
        Xt[[m]][t, , ] <- global_main + local_main + global_common + local_common + error_global + error_local + error_idiosyncratic
        if (store_components) {
          components$global_main[[m]][t, , ] <- global_main
          components$local_main[[m]][t, , ] <- local_main
          components$global_common[[m]][t, , ] <- global_common
          components$local_common[[m]][t, , ] <- local_common
          components$error_global[[m]][t, , ] <- error_global
          components$error_local[[m]][t, , ] <- error_local
          components$error_idiosyncratic[[m]][t, , ] <- error_idiosyncratic
        }
      }
      if (any(!is.finite(Xt[[m]]))) stop("Nonfinite generated observations in group ", if (is.null(group_names)) m else group_names[m], ".")
    }
    list(
      Xt = Xt, dimensions = dimensions, rank = .new_mmefm_rank(ranks),
      main_effect = list(
        direct = list(c = c, a = a, b = b),
        global = list(mu = mu, A1 = A1, alpha = alpha, B1 = B1, beta = beta),
        local = list(mu = mu_local, A2 = A2, alpha = alpha_local, B2 = B2, beta = beta_local)
      ),
      common_component = list(global = list(Q = Q, J = J, G = G, active = active), local = list(R = R, C = C, F = F)),
      error = list(global = list(Q = Q_e, J = J_e, G = G_e), local = list(R = R_e, C = C_e, F = F_e), sigma = sigma),
      memory = memory, innovation = list(family = innovation, df = df), components = components
    )
  })
}
