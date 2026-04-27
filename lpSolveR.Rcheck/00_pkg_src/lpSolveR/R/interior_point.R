#' Solve an LP Using the Primal-Dual Interior Point Method
#'
#' Implements a log-barrier interior point method with simple barrier
#' parameter reduction to solve a linear program. The problem is first
#' converted to standard form.
#'
#' @param obj Numeric vector of objective coefficients (length n).
#' @param A_ub Numeric matrix of <= inequality constraint coefficients (m x n).
#'   Default \code{NULL}.
#' @param b_ub Numeric vector of <= inequality RHS values (length m).
#'   Default \code{NULL}.
#' @param A_eq Numeric matrix of equality constraint coefficients (p x n).
#'   Default \code{NULL}.
#' @param b_eq Numeric vector of equality RHS values (length p).
#'   Default \code{NULL}.
#' @param maximize Logical. Maximise the objective? Default \code{FALSE}.
#' @param mu Initial barrier parameter. Default \code{10}.
#' @param tol Convergence tolerance for the duality gap. Default \code{1e-8}.
#' @param max_iter Maximum number of outer iterations. Default \code{200}.
#'
#' @return A named list:
#' \describe{
#'   \item{solution}{Numeric vector of optimal variable values.}
#'   \item{obj_value}{Optimal objective value.}
#'   \item{status}{One of \code{"optimal"} or \code{"max_iter_reached"}.}
#'   \item{iterations}{Number of iterations performed.}
#'   \item{runtime}{Elapsed time in seconds.}
#'   \item{duality_gap_history}{Numeric vector of duality gap at each iteration.}
#' }
#'
#' @examples
#' obj  <- c(2, 3)
#' A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' res  <- interior_point(obj, A_ub, b_ub)
#' res$obj_value
#'
#' @export
interior_point <- function(obj,
                            A_ub    = NULL, b_ub = NULL,
                            A_eq    = NULL, b_eq = NULL,
                            maximize = FALSE,
                            mu       = 10,
                            tol      = 1e-8,
                            max_iter = 200) {
  start_time <- proc.time()["elapsed"]

  sf <- lp_to_standard_form(obj, A_ub, b_ub, A_eq, b_eq,
                             maximize = maximize)
  A_mat  <- sf$A
  b_vec  <- sf$b
  c_vec  <- sf$c
  n_tot  <- ncol(A_mat)
  m      <- nrow(A_mat)
  n_orig <- sf$n_original

  # ---- Starting point -------------------------------------------------------
  # Use a simple interior point: x = e, compute s, y from KKT residuals.
  x <- rep(1, n_tot)

  # y via normal equations: A A^T y ≈ A c
  AAT <- A_mat %*% t(A_mat) + diag(1e-6, m)
  y   <- tryCatch(
    as.vector(solve(AAT, as.vector(A_mat %*% c_vec))),
    error = function(e) rep(0, m)
  )
  s   <- pmax(c_vec - as.vector(t(A_mat) %*% y), 1e-4)

  duality_gap_history <- numeric(0)
  status     <- "max_iter_reached"
  mu_current <- mu

  for (iter in seq_len(max_iter)) {
    # ---- Duality gap --------------------------------------------------------
    xs   <- x * s
    gap  <- sum(xs) / n_tot
    if (!is.finite(gap)) break   # numerical failure — stop cleanly
    duality_gap_history <- c(duality_gap_history, gap)
    if (gap < tol) { status <- "optimal"; break }

    # ---- Residuals ----------------------------------------------------------
    r_b <- b_vec - as.vector(A_mat %*% x)          # primal infeasibility
    r_c <- c_vec - as.vector(t(A_mat) %*% y) - s   # dual infeasibility

    # ---- Solve normal equations for dy  -------------------------------------
    # ( A diag(x/s) A^T ) dy = r_b + A diag(x/s)(r_c + mu/x)
    xs_ratio <- pmin(x / pmax(s, 1e-14), 1e12)     # clamp for stability
    rhs_vec  <- r_c + mu_current / pmax(x, 1e-14)
    M_mat    <- A_mat %*% (xs_ratio * t(A_mat)) + diag(1e-10, m)
    rhs_dy   <- r_b + as.vector(A_mat %*% (xs_ratio * rhs_vec))

    dy <- tryCatch(
      as.vector(solve(M_mat, rhs_dy)),
      error = function(e) rep(0, m)
    )
    ds <- r_c - as.vector(t(A_mat) %*% dy)
    dx <- (mu_current - xs - x * ds) / pmax(s, 1e-14)

    # Guard against non-finite directions
    if (!all(is.finite(dx)) || !all(is.finite(ds))) {
      dx[!is.finite(dx)] <- 0
      ds[!is.finite(ds)] <- 0
    }

    # ---- Step lengths -------------------------------------------------------
    alpha_p <- .safe_step(x, dx, 0.99)
    alpha_d <- .safe_step(s, ds, 0.99)

    x <- pmax(x + alpha_p * dx, 1e-14)
    y <- y + alpha_d * dy
    s <- pmax(s + alpha_d * ds, 1e-14)

    # ---- Reduce barrier parameter -------------------------------------------
    mu_current <- max(mu_current * 0.2, tol * 0.1)
  }

  solution <- x[1:n_orig]
  obj_val  <- sum(obj * solution)

  list(
    solution            = solution,
    obj_value           = obj_val,
    status              = status,
    iterations          = length(duality_gap_history),
    runtime             = as.numeric(proc.time()["elapsed"] - start_time),
    duality_gap_history = duality_gap_history
  )
}

.safe_step <- function(v, dv, fraction = 0.99) {
  neg_idx <- which(dv < 0 & is.finite(dv))
  if (length(neg_idx) == 0) return(1)
  min(fraction * min(-v[neg_idx] / dv[neg_idx]), 1)
}
