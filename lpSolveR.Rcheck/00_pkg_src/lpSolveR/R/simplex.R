#' Solve an LP Using the Simplex Method
#'
#' Solves a linear program using the full tableau simplex method with optional
#' Bland's rule to prevent cycling. The problem is first converted to standard
#' form via \code{\link{lp_to_standard_form}}.
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
#' @param tol Numeric tolerance for optimality and feasibility checks.
#'   Default \code{1e-8}.
#' @param max_iter Maximum number of pivot iterations. Default \code{1000}.
#' @param bland Logical. Use Bland's rule (smallest index entering variable) to
#'   prevent cycling? Default \code{TRUE}.
#'
#' @return A named list:
#' \describe{
#'   \item{solution}{Numeric vector of optimal variable values (original vars).}
#'   \item{obj_value}{Optimal objective value (in original sense).}
#'   \item{status}{One of \code{"optimal"}, \code{"infeasible"},
#'     \code{"unbounded"}, \code{"max_iter_reached"}.}
#'   \item{iterations}{Number of pivot steps performed.}
#'   \item{runtime}{Elapsed time in seconds.}
#'   \item{basis}{Integer vector of basic variable indices at termination.}
#'   \item{tableau_history}{List of objective values after each pivot.}
#' }
#'
#' @examples
#' obj  <- c(2, 3)
#' A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' res  <- simplex_method(obj, A_ub, b_ub)
#' res$obj_value
#' res$solution
#'
#' @export
simplex_method <- function(obj,
                            A_ub    = NULL, b_ub = NULL,
                            A_eq    = NULL, b_eq = NULL,
                            maximize = FALSE,
                            tol      = 1e-8,
                            max_iter = 1000,
                            bland    = TRUE) {
  start_time <- proc.time()["elapsed"]

  sf <- lp_to_standard_form(obj, A_ub, b_ub, A_eq, b_eq,
                             maximize = maximize)
  c_vec <- sf$c
  A     <- sf$A
  b     <- sf$b
  n_tot <- ncol(A)
  m     <- nrow(A)
  n_orig <- sf$n_original
  n_art  <- sf$n_artificial

  # Ensure b >= 0 (should be guaranteed by standard form builder)
  for (i in seq_len(m)) {
    if (b[i] < -tol) {
      A[i, ] <- -A[i, ]
      b[i]   <- -b[i]
    }
  }

  # Initial basis: prefer slack/artificial variables (last columns)
  # Slacks occupy columns n_orig+1 .. n_orig+n_slack+n_surplus
  # Artificials occupy the last n_art columns
  n_non_art <- n_tot - n_art
  basis <- integer(m)

  # Try to find identity columns for each row
  for (i in seq_len(m)) {
    found <- FALSE
    # Check artificials first (rightmost) - only if they exist
    if (n_art > 0) {
      for (j in seq(n_non_art + 1, n_tot)) {
        col <- A[, j]
        if (abs(col[i] - 1) < tol && sum(abs(col)) - 1 < tol) {
          basis[i] <- j
          found <- TRUE
          break
        }
      }
    }
    if (!found && n_non_art > n_orig) {
      # Try slacks / surplus
      for (j in seq(n_orig + 1, n_non_art)) {
        col <- A[, j]
        if (abs(col[i] - 1) < tol && sum(abs(col)) - 1 < tol) {
          basis[i] <- j
          found <- TRUE
          break
        }
      }
    }
    if (!found) basis[i] <- i  # fallback
  }

  # Build initial tableau: [A | b] with reduced costs
  # Tableau rows: constraint rows + objective row
  # T[1:m, ] = A augmented with b; T[m+1, ] = reduced cost row
  tab <- cbind(A, b)  # m x (n_tot+1)
  obj_row <- c(c_vec, 0)  # will be updated after pivoting in basis

  # Re-express objective in terms of non-basic variables
  for (i in seq_len(m)) {
    bv <- basis[i]
    obj_row <- obj_row - obj_row[bv] * tab[i, ]
  }

  iter <- 0
  tableau_history <- numeric(0)

  repeat {
    iter <- iter + 1
    if (iter > max_iter) {
      status <- "max_iter_reached"
      break
    }

    reduced_costs <- obj_row[1:n_tot]

    # Check optimality
    if (bland) {
      entering <- which(reduced_costs < -tol)[1]
    } else {
      min_rc <- min(reduced_costs)
      if (min_rc >= -tol) { status <- "optimal"; break }
      entering <- which.min(reduced_costs)
    }
    if (is.na(entering)) { status <- "optimal"; break }

    # Min-ratio test (Bland: smallest ratio, ties broken by smallest index)
    col_enter <- tab[, entering]
    ratios <- ifelse(col_enter > tol, tab[, n_tot + 1] / col_enter, Inf)

    if (all(ratios == Inf)) {
      status <- "unbounded"
      break
    }

    min_ratio <- min(ratios)
    leaving_rows <- which(abs(ratios - min_ratio) < tol)
    if (bland && length(leaving_rows) > 1) {
      # Bland: leave row whose basic variable has smallest index
      leaving <- leaving_rows[which.min(basis[leaving_rows])]
    } else {
      leaving <- leaving_rows[1]
    }

    # Pivot
    pivot_val <- tab[leaving, entering]
    tab[leaving, ] <- tab[leaving, ] / pivot_val

    for (i in seq_len(m)) {
      if (i != leaving) {
        tab[i, ] <- tab[i, ] - tab[i, entering] * tab[leaving, ]
      }
    }
    obj_row <- obj_row - obj_row[entering] * tab[leaving, ]

    basis[leaving] <- entering
    tableau_history <- c(tableau_history, -obj_row[n_tot + 1])
  }

  if (iter <= max_iter && !exists("status")) status <- "optimal"

  # Check if artificials are still in basis (infeasibility)
  if (status == "optimal" && n_art > 0) {
    art_in_basis <- basis > n_non_art
    if (any(art_in_basis)) {
      # Check if their values are non-zero
      art_vals <- tab[art_in_basis, n_tot + 1]
      if (any(abs(art_vals) > tol)) {
        status <- "infeasible"
      }
    }
  }

  # Extract solution
  x_full <- numeric(n_tot)
  for (i in seq_len(m)) x_full[basis[i]] <- tab[i, n_tot + 1]

  solution <- x_full[1:n_orig]

  obj_val <- if (status == "optimal") {
    val <- sum(obj * solution)
    val
  } else NA_real_

  runtime <- proc.time()["elapsed"] - start_time

  list(
    solution        = solution,
    obj_value       = obj_val,
    status          = status,
    iterations      = iter - (status != "max_iter_reached"),
    runtime         = as.numeric(runtime),
    basis           = basis,
    tableau_history = tableau_history
  )
}
