.standardize_detection_residuals <- function(check_Y) {
  scale <- numeric(length(check_Y))
  names(scale) <- names(check_Y)
  check_Z <- check_Y
  for (m in seq_along(check_Y)) {
    scale[m] <- sqrt(mean(check_Y[[m]]^2))
    if (!is.finite(scale[m]) || scale[m] <= 0) {
      group <- if (is.null(names(check_Y))) m else names(check_Y)[m]
      stop("Residual RMS scale must be finite and strictly positive in group ", group, ".")
    }
    check_Z[[m]] <- check_Y[[m]] / scale[m]
  }
  list(check_Z = check_Z, scale = scale)
}

.global_detection_statistic <- function(check_Z, rank, Q_hat, J_hat) {
  M <- length(check_Z)
  TT <- dim(check_Z[[1L]])[1L]
  Z <- vector("list", M)
  denominator <- numeric(M)
  for (m in seq_len(M)) {
    p <- dim(check_Z[[m]])[2L]
    q <- dim(check_Z[[m]])[3L]
    row_cov <- matrix(0, p, p)
    col_cov <- matrix(0, q, q)
    Z[[m]] <- matrix(0, TT, as.numeric(rank$kr) * rank$kc)
    for (t in seq_len(TT)) {
      Yt <- matrix(check_Z[[m]][t, , , drop = FALSE], p, q)
      Z[[m]][t, ] <- as.vector(crossprod(Q_hat[[m]], Yt) %*% J_hat[[m]])
      row_cov <- row_cov + tcrossprod(Yt) / TT
      col_cov <- col_cov + crossprod(Yt) / TT
    }
    # The denominator uses total common spaces, without a global/local split.
    U_hat <- .leading_eigenvectors(row_cov, rank$kr + rank$kr_m[m])
    V_hat <- .leading_eigenvectors(col_cov, rank$kc + rank$kc_m[m])
    energy <- matrix(0, ncol(U_hat), ncol(V_hat))
    for (t in seq_len(TT)) {
      Yt <- matrix(check_Z[[m]][t, , , drop = FALSE], p, q)
      energy <- energy + (crossprod(U_hat, Yt) %*% V_hat)^2 / TT
    }
    denominator[m] <- max(energy)
    if (!is.finite(denominator[m]) || denominator[m] <= 0) {
      group <- if (is.null(names(check_Z))) m else names(check_Z)[m]
      stop("Detection denominator must be finite and strictly positive in group ", group, ".")
    }
  }
  SG <- numeric(M)
  names(SG) <- names(check_Z)
  for (m in seq_len(M)) {
    numerator <- 0
    for (n in seq_len(M)) {
      if (n != m) {
        numerator <- numerator + max((crossprod(Z[[m]], Z[[n]]) / TT)^2)
      }
    }
    SG[m] <- numerator / denominator[m]^2
    if (!is.finite(SG[m])) {
      stop("Nonfinite detection statistic in group ", if (is.null(names(check_Z))) m else names(check_Z)[m], ".")
    }
  }
  list(statistic = SG, second_largest = unname(sort(SG, decreasing = TRUE)[2L]))
}

.circular_shift_detection_groups <- function(check_Z, reference_group, shifts = NULL) {
  TT <- dim(check_Z[[1L]])[1L]
  shifted <- check_Z
  for (m in seq_along(check_Z)) {
    if (m != reference_group) {
      shift <- if (is.null(shifts)) sample.int(TT, 1L) - 1L else shifts[m]
      index <- 1L + ((seq_len(TT) - 1L + as.numeric(shift)) %% TT)
      shifted[[m]] <- check_Z[[m]][index, , , drop = FALSE]
    }
  }
  shifted
}

.bootstrap_global_detection_one <- function(job, inputs = get(".mmefm_detection_inputs", envir = .GlobalEnv),
  verbose = FALSE) {
  messages <- character()
  precomputed <- !is.null(job$initial_candidate_indices)
  draw <- tryCatch(withCallingHandlers({
    shifted <- .circular_shift_detection_groups(inputs$check_Z, inputs$reference_group, job$shifts)
    loadings <- .initial_common_loadings(shifted, inputs$rank$kr, inputs$rank$kc, inputs$K0,
                                         inputs$max_iter, inputs$tol, verbose, job$initial_candidate_indices)
    statistic <- .global_detection_statistic(shifted, inputs$rank, loadings$Q_hat, loadings$J_hat)
    list(second_largest = statistic$second_largest, converged = loadings$converged)
  }, warning = function(w) {
    if (precomputed) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  }), error = function(e) {
    if (!precomputed) {
      stop("Bootstrap replication ", job$b, " failed: ", conditionMessage(e), call. = FALSE)
    }
    list(error = conditionMessage(e))
  })
  draw$warnings <- messages
  draw
}

.bootstrap_global_detection <- function(check_Z, rank, B, reference_group, K0, max_iter, tol, verbose,
                                         parallel = FALSE, num.cores = 4L) {
  inputs <- list(check_Z = check_Z, rank = rank, reference_group = reference_group,
                 K0 = K0, max_iter = max_iter, tol = tol)
  second_largest <- numeric(B)
  nonconverged <- 0L
  if (parallel && num.cores > 1L && B > 1L) {
    M <- length(check_Z)
    TT <- dim(check_Z[[1L]])[1L]
    # Generate shifts and half-panels on the master in serial RNG order.
    jobs <- vector("list", B)
    for (b in seq_len(B)) {
      shifts <- integer(M)
      for (m in seq_len(M)) {
        if (m != reference_group) {
          shifts[m] <- sample.int(TT, 1L) - 1L
        }
      }
      indices <- vector("list", M)
      for (m in seq_len(M)) {
        q_m <- dim(check_Z[[m]])[3L]
        indices[[m]] <- vector("list", K0)
        for (k in seq_len(K0)) {
          indices[[m]][[k]] <- sample(q_m, floor(q_m / 2))
        }
      }
      jobs[[b]] <- list(b = b, shifts = shifts, initial_candidate_indices = indices)
    }
    package_path <- getNamespaceInfo(asNamespace("MMEFM"), "path")
    if (!file.exists(file.path(package_path, "Meta", "package.rds"))) {
      stop("PSOCK bootstrap requires the current installed MMEFM package; install this source before parallel development testing.")
    }
    cl <- parallel::makePSOCKcluster(min(B, num.cores))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    worker_setup <- function(package_path, library_paths) {
      .libPaths(c(dirname(package_path), library_paths))
      namespace <- loadNamespace("MMEFM", lib.loc = dirname(package_path))
      if (normalizePath(getNamespaceInfo(namespace, "path")) != normalizePath(package_path)) {
        stop("PSOCK worker loaded a different MMEFM installation.")
      }
      NULL
    }
    # Resolve the installation before deserializing package functions; do not ship the master's call frame.
    environment(worker_setup) <- baseenv()
    parallel::clusterCall(cl, worker_setup, package_path, .libPaths())
    # Fixed arrays and controls are sent once per worker, separately from compact jobs.
    .mmefm_detection_inputs <- inputs
    parallel::clusterExport(cl, ".mmefm_detection_inputs", envir = environment())
    if (verbose) {
      cat("Running", B, "bootstrap replications on", min(B, num.cores), "PSOCK workers.\n")
    }
    draws <- parallel::parLapplyLB(cl, jobs, .bootstrap_global_detection_one)
    for (b in seq_len(B)) {
      for (message in draws[[b]]$warnings) {
        warning("Bootstrap replication ", b, ": ", message, call. = FALSE)
      }
    }
    for (b in seq_len(B)) {
      draw <- draws[[b]]
      if (!is.null(draw$error)) {
        stop("Bootstrap replication ", b, " failed: ", draw$error, call. = FALSE)
      }
      second_largest[b] <- draw$second_largest
      if (!draw$converged) {
        nonconverged <- nonconverged + 1L
      }
    }
    if (verbose) {
      cat("Parallel bootstrap complete.\n")
    }
  } else {
    for (b in seq_len(B)) {
      draw <- .bootstrap_global_detection_one(list(b = b), inputs, if (parallel) FALSE else verbose)
      second_largest[b] <- draw$second_largest
      if (!draw$converged) {
        nonconverged <- nonconverged + 1L
      }
    }
  }
  list(second_largest = second_largest, nonconverged = nonconverged)
}

#' Detect and screen global common factors in MMEFM
#'
#' Test whether a global common component is shared by at least two groups, and screen the groups that share it. Group statistics are calibrated by a circular-shift bootstrap.
#' @param Xt A list of finite `T x p_m x q_m` arrays as in [est_MMEFM()], or an `mmefm_fit`. Arrays are reduced to direct additive residuals `check_Y`; a fit supplies its stored `check_Y`.
#' @param rank NULL to select ranks for array input or use a fit's ranks. Alternatively, supply a complete eight-field rank list or `mmefm_rank` as in [est_MMEFM()]. Supplied ranks override fit ranks. Rank values are validated, and any rank-object diagnostics are retained.
#' @param B Positive whole-number bootstrap size. The default is 199. Values such as 3 are suitable only for illustrative examples, not reliable calibration.
#' @param alpha Finite significance level strictly between zero and one, used for the bootstrap cutoff.
#' @param reference_group Integer index from 1 through M identifying the group held fixed during bootstrap shifts.
#' @param parallel One nonmissing logical. FALSE runs the bootstrap serially. TRUE uses PSOCK workers on Windows, macOS, and Linux when B and `num.cores` both exceed one.
#' @param num.cores Positive whole-number maximum worker count. At most `min(B, num.cores)` workers are used; one worker or B = 1 gives serial execution. Stay within the available resource allocation.
#' @inheritParams select_MMEFM_rank
#' @param seed Nonnegative whole number within R's integer range controlling projection directions and bootstrap shifts. The caller's RNG state is restored on return or error.
#' @details
#' Detection concerns global common matrix-factor interactions, not global main effects. Each group's direct residual is divided by its RMS scale without further centering. Initial global loading spaces are re-estimated from these standardized residuals, including for fit input.
#'
#' The group statistic \eqn{S_{G,m}} measures cross-group dependence relative to within-group variation. Its numerator sums the largest squared cross-group moments of projected global-core coordinates over other groups. Its denominator is the square of the largest within-group projected-coordinate mean square, using total global plus local common loading dimensions. These moments are not additionally demeaned over time.
#'
#' The existence test uses the second-largest group statistic because a global common component must involve at least two groups. The bootstrap holds the reference group fixed and independently circularly shifts each other group's time series by a uniformly drawn offset from 0 through T-1. Loading spaces and statistics are re-estimated at fixed ranks for each replication. This retains within-group temporal dependence while disrupting synchronous cross-group alignment.
#'
#' One Type-8 empirical `(1 - alpha)` quantile of bootstrap second-largest statistics supplies the common cutoff. Existence requires the observed second-largest statistic to strictly exceed it. Final screening retains groups strictly above the same cutoff only if existence is detected and at least two groups qualify; otherwise the selected set is empty.
#'
#' Whole-sample circular shifts assume a single stationary segment. They are not justified across arbitrary structural breaks or piecewise nonstationary series. The model's allowance for nonstationary main effects does not remove this bootstrap restriction.
#'
#' Invalid residual scales, denominators, and numerical failures error without repair or discarded bootstrap draws. Bootstrap errors identify the replication. Projection non-convergence retains final iterates and warns for the observed statistic or, in aggregate, for bootstrap replications. The caller's RNG state is restored on success or error; automatic rank selection uses a separate scoped seed.
#'
#' Parallel execution requires the current installed MMEFM namespace. Serial and parallel calls use matching replication-specific randomness under the same seed; workers are cleaned up on completion or error. Process startup and memory overhead can make small bootstraps slower. Avoid nested worker pools when replications are already parallelized externally.
#' @return An `mmefm_global_detection` object with `call`, group `statistic`, observed `second_largest`, `cutoff`, and logical `exists`. `initial_selected` marks groups exceeding the cutoff; `selected` applies the existence and at-least-two-groups requirements. Group vectors preserve optional names.
#'
#' `scale` contains group residual RMS values, and `rank` contains the validated `mmefm_rank`. `bootstrap` contains `B`, `alpha`, `reference_group`, and the length-B vector `second_largest` used for calibration. Loading estimates and residual arrays are not returned.
#' @seealso [est_MMEFM()], [select_MMEFM_rank()], [gen_MMEFM()]
#' @examples
#' rank <- list(r1 = 1L, l1 = 1L, r2 = rep(1L, 3), l2 = rep(1L, 3),
#'              kr = 1L, kc = 1L, kr_m = rep(1L, 3), kc_m = rep(1L, 3))
#' simulation <- gen_MMEFM(T = 40L, p = c(a = 6L, b = 7L, c = 8L), q = c(a = 7L, b = 8L, c = 6L),
#'                         rank = rank, burn = 20L, filter_length = 40L, seed = 2026L)
#' # B = 3 illustrates the call; it does not provide reliable calibration.
#' detection <- detect_MMEFM_global(simulation$Xt, rank = simulation$rank, B = 3L,
#'                                  K0 = 8L, max_iter = 50L, parallel = FALSE, seed = 2026L)
#' detection
#' detection$selected
#' @export
detect_MMEFM_global <- function(Xt, rank = NULL, B = 199L, alpha = 0.05, reference_group = 1L,
                               K0 = 20L, max_iter = 20L, tol = 1e-4, seed = 2026L, verbose = FALSE,
                               parallel = FALSE, num.cores = 4L) {
  call <- match.call()
  fit_input <- inherits(Xt, "mmefm_fit")
  data_info <- .validate_Xt(if (fit_input) Xt$check_Y else Xt)
  controls <- list(B = B, reference_group = reference_group, K0 = K0, max_iter = max_iter, seed = seed,
    num.cores = num.cores)
  for (name in names(controls)) {
    value <- controls[[name]]
    minimum <- if (name == "seed") 0L else 1L
    maximum <- if (name == "reference_group") min(data_info$M, .Machine$integer.max) else .Machine$integer.max
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != floor(value) ||
      value < minimum || value > maximum) {
      stop(name, " must be one finite whole number between ", minimum, " and ", maximum, ".")
    }
  }
  B <- as.integer(B)
  reference_group <- as.integer(reference_group)
  K0 <- as.integer(K0)
  max_iter <- as.integer(max_iter)
  seed <- as.integer(seed)
  num.cores <- as.integer(num.cores)
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be one nonmissing logical.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite number strictly between zero and one.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one nonmissing logical.")
  }
  if (fit_input) {
    check_Y <- Xt$check_Y
    if (is.null(rank)) {
      rank <- Xt$rank
    }
  } else {
    check_Y <- .main_effect_moments(Xt)$check_Y
    if (is.null(rank)) {
      rank <- select_MMEFM_rank(Xt, K0, max_iter, tol, seed, verbose)
    }
  }
  ranks <- .validate_rank(rank, data_info)
  rank <- .new_mmefm_rank(ranks, if (inherits(rank, "mmefm_rank")) rank$diagnostics else NULL)

  # RMS scaling retains the direct main-effect residual geometry.
  standardized <- .standardize_detection_residuals(check_Y)
  computed <- .with_preserved_seed(seed, {
    loadings <- .initial_common_loadings(standardized$check_Z, ranks$kr, ranks$kc, K0, max_iter, tol, verbose)
    if (!loadings$converged) {
      warning("Observed detection projection reached max_iter without convergence.", call. = FALSE)
    }
    observed <- .global_detection_statistic(standardized$check_Z, ranks, loadings$Q_hat, loadings$J_hat)
    bootstrap <- .bootstrap_global_detection(standardized$check_Z, ranks, B, reference_group, K0, max_iter,
      tol, verbose, parallel, num.cores)
    if (bootstrap$nonconverged > 0L) {
      warning(bootstrap$nonconverged, " of ", B,
        " bootstrap projections reached max_iter without convergence.", call. = FALSE)
    }
    list(observed = observed, bootstrap = bootstrap)
  })

  # One Type-8 cutoff calibrates both existence and group screening.
  cutoff <- stats::quantile(computed$bootstrap$second_largest, 1 - alpha, names = FALSE, type = 8L)
  observed <- computed$observed
  exists <- observed$second_largest > cutoff
  initial_selected <- observed$statistic > cutoff
  selected <- initial_selected
  if (!exists || sum(initial_selected) < 2L) {
    selected[] <- FALSE
  }
  structure(
    list(
      call = call, statistic = observed$statistic, second_largest = observed$second_largest,
      cutoff = cutoff, exists = exists, initial_selected = initial_selected, selected = selected,
      scale = standardized$scale, rank = rank,
      bootstrap = list(
        B = B, alpha = alpha, reference_group = reference_group,
        second_largest = computed$bootstrap$second_largest
      )
    ),
    class = "mmefm_global_detection"
  )
}

#' Print global common-factor detection results
#'
#' Display the second-largest statistic, bootstrap cutoff, calibration level, existence decision, and selected groups.
#' @param x An `mmefm_global_detection` returned by [detect_MMEFM_global()].
#' @param ... Further arguments, currently unused.
#' @return The input object invisibly.
#' @seealso [detect_MMEFM_global()]
#' @export
print.mmefm_global_detection <- function(x, ...) {
  cat(sprintf("MMEFM global detection: %d groups\n", length(x$statistic)))
  cat(sprintf("Second-largest S_G: %.6g; cutoff: %.6g; alpha: %.3g\n", x$second_largest, x$cutoff, x$bootstrap$alpha))
  cat("Global common factors detected:", if (x$exists) "yes" else "no", "\n")
  groups <- which(x$selected)
  labels <- if (is.null(names(x$selected))) as.character(groups) else names(x$selected)[groups]
  cat("Selected groups:", if (length(labels)) paste(labels, collapse = ", ") else "none", "\n")
  invisible(x)
}
