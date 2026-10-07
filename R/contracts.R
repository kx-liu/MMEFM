# Canonical axes are time, rows, columns; validation leaves the arrays untouched.
.validate_Xt <- function(Xt) {
  if (!is.list(Xt) || is.data.frame(Xt)) stop("Xt must be a list of group arrays.")
  M <- length(Xt)
  if (M < 2L) stop("Xt must contain at least two groups.")
  group_names <- names(Xt)
  if (!is.null(group_names) &&
      (anyNA(group_names) || any(!nzchar(group_names)) || anyDuplicated(group_names))) {
    stop("Xt group names must be nonempty, nonmissing, and unique.")
  }

  TT <- NULL
  p <- q <- integer(M)
  for (m in seq_len(M)) {
    x <- Xt[[m]]
    dims <- dim(x)
    if (!is.numeric(x) || length(dims) != 3L) {
      stop(sprintf("Xt group %d must be a numeric three-dimensional T x p_m x q_m array.", m))
    }
    if (dims[1L] < 1L || dims[2L] < 2L || dims[3L] < 2L) {
      stop(sprintf("Xt group %d must have T >= 1, p_m >= 2, and q_m >= 2.", m))
    }
    if (is.null(TT)) TT <- dims[1L]
    if (dims[1L] != TT) stop("All Xt groups must have the same time dimension T.")
    if (any(!is.finite(x))) stop(sprintf("Xt group %d must contain only finite entries.", m))
    p[m] <- dims[2L]
    q[m] <- dims[3L]
  }
  names(p) <- names(q) <- group_names
  list(M = M, T = TT, p = p, q = q, group_names = group_names)
}

# data_info is the metadata returned by .validate_Xt(); ranks are complete.
.validate_rank <- function(rank, data_info) {
  components <- c("r1", "l1", "r2", "l2", "kr", "kc", "kr_m", "kc_m")
  if (!is.list(rank) || is.data.frame(rank) || length(rank) != length(components) ||
      is.null(names(rank)) || anyNA(names(rank)) || anyDuplicated(names(rank)) ||
      !setequal(names(rank), components)) {
    stop("rank must be a complete named list containing exactly r1, l1, r2, l2, kr, kc, kr_m, and kc_m.")
  }

  ranks <- rank[components]
  for (component in components) {
    values <- ranks[[component]]
    global <- component %in% c("r1", "l1", "kr", "kc")
    required_length <- if (global) 1L else data_info$M
    if (!is.numeric(values) || !is.null(dim(values)) || length(values) != required_length ||
        any(!is.finite(values)) || any(values <= 0) || any(values != floor(values))) {
      stop(sprintf("rank$%s must be a vector of %d positive whole-number value(s).", component, required_length))
    }

    if (!global) {
      value_names <- names(values)
      if (!is.null(value_names)) {
        if (anyNA(value_names) || any(!nzchar(value_names)) || anyDuplicated(value_names)) {
          stop(sprintf("rank$%s group names must be nonempty, nonmissing, and unique.", component))
        }
        if (!is.null(data_info$group_names)) {
          if (!setequal(value_names, data_info$group_names)) {
            stop(sprintf("rank$%s names must match the Xt group names exactly.", component))
          }
          values <- values[match(data_info$group_names, value_names)]
        }
      }
    }

    row_side <- component %in% c("r1", "r2", "kr", "kr_m")
    dimensions <- if (row_side) data_info$p else data_info$q
    limit <- if (global) min(dimensions) else dimensions
    if (any(values >= limit)) {
      requirement <- if (global) "the minimum group dimension" else "each corresponding group dimension"
      stop(sprintf("rank$%s must be smaller than %s on its %s side.",
                   component, requirement, if (row_side) "row" else "column"))
    }
    values <- as.integer(values)
    if (!global) names(values) <- data_info$group_names
    ranks[[component]] <- values
  }
  ranks
}

# The constructor receives canonical ranks and does not repeat validation.
.new_mmefm_rank <- function(ranks, diagnostics = NULL) {
  structure(list(r1 = ranks$r1, l1 = ranks$l1, r2 = ranks$r2, l2 = ranks$l2,
                 kr = ranks$kr, kc = ranks$kc, kr_m = ranks$kr_m, kc_m = ranks$kc_m,
                 diagnostics = diagnostics), class = "mmefm_rank")
}
