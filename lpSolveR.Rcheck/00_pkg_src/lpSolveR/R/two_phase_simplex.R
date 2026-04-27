#' Solve an LP Using the Two-Phase Simplex Method
#'
#' Implements the two-phase simplex method. Phase 1 minimises the sum of
#' artificial variables to find a basic feasible solution (BFS). Phase 2 then
#' optimises the original objective from that BFS.
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
#' @param tol Numeric tolerance. Default \code{1e-8}.
#' @param max_iter Maximum pivot iterations per phase. Default \code{1000}.
#'
#' @return A named list:
#' \describe{
#'   \item{solution}{Numeric vector of optimal variable values.}
#'   \item{obj_value}{Optimal objective value.}
#'   \item{status}{One of \code{"optimal"}, \code{"infeasible"},
#'     \code{"unbounded"}, \code{"max_iter_reached"}.}
#'   \item{iterations}{Total pivot steps (both phases).}
#'   \item{runtime}{Elapsed time in seconds.}
#'   \item{basis}{Basic variable indices at termination.}
#'   \item{tableau_history}{Objective values after each Phase 2 pivot.}
#'   \item{phase1_iterations}{Pivot steps in Phase 1.}
#'   \item{phase2_iterations}{Pivot steps in Phase 2.}
#' }
#'
#' @examples
#' obj  <- c(1, 2)
#' A_eq <- matrix(c(1, 1), nrow = 1)
#' b_eq <- 3
#' res  <- two_phase_simplex(obj, A_eq = A_eq, b_eq = b_eq)
#' res$obj_value
#'
#' @export
two_phase_simplex <- function(obj,
                               A_ub    = NULL, b_ub = NULL,
                               A_eq    = NULL, b_eq = NULL,
                               maximize = FALSE,
                               tol      = 1e-8,
                               max_iter = 1000) {
  start_time <- proc.time()["elapsed"]

  sf <- lp_to_standard_form(obj, A_ub, b_ub, A_eq, b_eq,
                             maximize = maximize)

  A     <- sf$A
  b     <- sf$b
  c_orig <- sf$c
  n_tot  <- ncol(A)
  m      <- nrow(A)
  n_orig <- sf$n_original
  n_art  <- sf$n_artificial
  n_non_art <- n_tot - n_art

  # Ensure b >= 0
  for (i in seq_len(m)) {
    if (b[i] < -tol) {
      A[i, ] <- -A[i, ]
      b[i]   <- -b[i]
    }
  }

  # ------------------------------------------------------------------ Phase 1
  # Objective: minimise sum of artificial variables
  c_p1 <- c(rep(0, n_non_art), rep(1, n_art))

  # Initial basis: artificials (last n_art columns) + slacks where available
  basis <- .init_basis(A, b, n_orig, n_non_art, n_tot, m, tol)

  tab <- cbind(A, b)
  obj_row <- c(c_p1, 0)
  for (i in seq_len(m)) {
    bv <- basis[i]
    obj_row <- obj_row - obj_row[bv] * tab[i, ]
  }

  iter_p1 <- 0
  status   <- "running"

  repeat {
    iter_p1 <- iter_p1 + 1
    if (iter_p1 > max_iter) { status <- "max_iter_reached"; break }

    rc <- obj_row[1:n_tot]
    entering <- which(rc < -tol)[1]
    if (is.na(entering)) { status <- "p1_done"; break }

    col_enter <- tab[, entering]
    ratios <- ifelse(col_enter > tol, tab[, n_tot + 1] / col_enter, Inf)
    if (all(ratios == Inf)) { status <- "unbounded"; break }

    min_r <- min(ratios)
    leaving_rows <- which(abs(ratios - min_r) < tol)
    leaving <- leaving_rows[which.min(basis[leaving_rows])]

    tab <- .pivot(tab, leaving, entering, m)
    obj_row <- obj_row - obj_row[entering] * tab[leaving, ]
    basis[leaving] <- entering
  }

  if (status == "max_iter_reached") {
    return(.make_result(NA, NA, "max_iter_reached", iter_p1, 0,
                        basis, numeric(0), iter_p1, 0, start_time))
  }

  # Check Phase 1 objective value (should be 0 for feasibility)
  p1_obj_val <- -obj_row[n_tot + 1]
  if (p1_obj_val > tol) {
    return(.make_result(rep(NA, n_orig), NA, "infeasible", iter_p1, 0,
                        basis, numeric(0), iter_p1, 0, start_time))
  }

  # Drive artificials out of basis if possible
  for (i in seq_len(m)) {
    if (basis[i] > n_non_art) {
      for (j in seq_len(n_non_art)) {
        if (abs(tab[i, j]) > tol) {
          tab  <- .pivot(tab, i, j, m)
          obj_row <- obj_row - obj_row[j] * tab[i, ]
          basis[i] <- j
          break
        }
      }
    }
  }

  # Remove artificial columns from tableau
  keep_cols <- c(seq_len(n_non_art), n_tot + 1)
  tab2   <- tab[, keep_cols, drop = FALSE]
  n_tot2 <- n_non_art

  # ------------------------------------------------------------------ Phase 2
  c_p2    <- c_orig[1:n_non_art]
  obj_row2 <- c(c_p2, 0)
  for (i in seq_len(m)) {
    bv <- basis[i]
    if (bv <= n_non_art) {
      obj_row2 <- obj_row2 - obj_row2[bv] * tab2[i, ]
    }
  }

  iter_p2 <- 0
  tableau_history <- numeric(0)
  status2 <- "running"

  repeat {
    iter_p2 <- iter_p2 + 1
    if (iter_p2 > max_iter) { status2 <- "max_iter_reached"; break }

    rc2 <- obj_row2[1:n_tot2]
    entering2 <- which(rc2 < -tol)[1]
    if (is.na(entering2)) { status2 <- "optimal"; break }

    col2 <- tab2[, entering2]
    ratios2 <- ifelse(col2 > tol, tab2[, n_tot2 + 1] / col2, Inf)
    if (all(ratios2 == Inf)) { status2 <- "unbounded"; break }

    min_r2 <- min(ratios2)
    leaving_rows2 <- which(abs(ratios2 - min_r2) < tol)
    leaving2 <- leaving_rows2[which.min(basis[leaving_rows2])]

    tab2 <- .pivot(tab2, leaving2, entering2, m)
    obj_row2 <- obj_row2 - obj_row2[entering2] * tab2[leaving2, ]
    basis[leaving2] <- entering2
    tableau_history <- c(tableau_history, -obj_row2[n_tot2 + 1])
  }

  if (status2 == "running") status2 <- "optimal"

  # Extract solution
  x_full <- numeric(n_non_art)
  for (i in seq_len(m)) {
    bv <- basis[i]
    if (bv <= n_non_art) x_full[bv] <- tab2[i, n_tot2 + 1]
  }

  solution  <- x_full[1:n_orig]
  obj_val   <- if (status2 == "optimal") sum(obj * solution) else NA_real_

  .make_result(solution, obj_val, status2,
               iter_p1 + iter_p2 - 1, iter_p2 - 1,
               basis, tableau_history,
               iter_p1 - 1, iter_p2 - 1, start_time)
}

# ---- helpers -----------------------------------------------------------------

.pivot <- function(tab, pivot_row, pivot_col, m) {
  tab[pivot_row, ] <- tab[pivot_row, ] / tab[pivot_row, pivot_col]
  for (i in seq_len(m)) {
    if (i != pivot_row) {
      tab[i, ] <- tab[i, ] - tab[i, pivot_col] * tab[pivot_row, ]
    }
  }
  tab
}

.init_basis <- function(A, b, n_orig, n_non_art, n_tot, m, tol) {
  basis <- integer(m)
  n_art <- n_tot - n_non_art
  for (i in seq_len(m)) {
    found <- FALSE
    # Prefer artificial variable columns (guaranteed identity)
    if (n_art > 0) {
      for (j in seq(n_non_art + 1, n_tot)) {
        col <- A[, j]
        if (abs(col[i] - 1) < tol && sum(abs(col)) - 1 < tol) {
          basis[i] <- j; found <- TRUE; break
        }
      }
    }
    # Fall back to slack columns
    if (!found && n_non_art > n_orig) {
      for (j in seq(n_orig + 1, n_non_art)) {
        col <- A[, j]
        if (abs(col[i] - 1) < tol && sum(abs(col)) - 1 < tol) {
          basis[i] <- j; found <- TRUE; break
        }
      }
    }
    if (!found) basis[i] <- if (n_art > 0) n_non_art + i else i
  }
  basis
}

.make_result <- function(solution, obj_value, status,
                          iterations, iter_p2_val,
                          basis, tableau_history,
                          phase1_iterations, phase2_iterations,
                          start_time) {
  list(
    solution          = solution,
    obj_value         = obj_value,
    status            = status,
    iterations        = max(0, iterations),
    runtime           = as.numeric(proc.time()["elapsed"] - start_time),
    basis             = basis,
    tableau_history   = tableau_history,
    phase1_iterations = max(0, phase1_iterations),
    phase2_iterations = max(0, phase2_iterations)
  )
}
