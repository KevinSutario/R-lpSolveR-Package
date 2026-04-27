#' Convert a General LP to Standard Form
#'
#' Converts a general linear program to standard form: minimize c^T x subject
#' to Ax = b, x >= 0. Slack variables are added for <= inequalities, surplus
#' and artificial variables for >= inequalities, and artificial variables for
#' equality constraints that lack an obvious basic feasible solution.
#'
#' @param obj Numeric vector of objective coefficients (length n).
#' @param A_ub Numeric matrix of <= inequality constraint coefficients (m1 x n).
#'   Default \code{NULL}.
#' @param b_ub Numeric vector of <= inequality RHS values (length m1).
#'   Default \code{NULL}.
#' @param A_eq Numeric matrix of equality constraint coefficients (m2 x n).
#'   Default \code{NULL}.
#' @param b_eq Numeric vector of equality RHS values (length m2).
#'   Default \code{NULL}.
#' @param lb Numeric vector of lower bounds for variables (length n).
#'   Default \code{NULL} (treated as 0).
#' @param ub Numeric vector of upper bounds for variables (length n).
#'   Default \code{NULL} (treated as Inf).
#' @param maximize Logical. If \code{TRUE}, the objective is negated to convert
#'   a maximisation problem to minimisation. Default \code{FALSE}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{c}{Objective coefficient vector in standard form.}
#'   \item{A}{Constraint matrix in standard form.}
#'   \item{b}{RHS vector in standard form (all non-negative after variable shift).}
#'   \item{n_original}{Number of original decision variables.}
#'   \item{n_slack}{Number of slack/surplus variables added.}
#'   \item{n_artificial}{Number of artificial variables added.}
#'   \item{variable_names}{Character vector of variable names.}
#'   \item{maximize}{Logical indicating if original problem was maximisation.}
#' }
#'
#' @examples
#' obj   <- c(2, 3)
#' A_ub  <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub  <- c(4, 6)
#' sf    <- lp_to_standard_form(obj, A_ub, b_ub)
#' sf$A
#' sf$variable_names
#'
#' @export
lp_to_standard_form <- function(obj,
                                 A_ub = NULL, b_ub = NULL,
                                 A_eq = NULL, b_eq = NULL,
                                 lb   = NULL, ub   = NULL,
                                 maximize = FALSE) {
  n <- length(obj)

  # Default bounds
  if (is.null(lb)) lb <- rep(0, n)
  if (is.null(ub)) ub <- rep(Inf, n)

  # Flip objective for maximisation
  c_vec <- if (maximize) -obj else obj

  # Shift variables so lower bounds are zero: x' = x - lb
  # Adjust RHS accordingly
  shift <- lb

  rows_ub  <- if (!is.null(A_ub)) nrow(A_ub) else 0
  rows_eq  <- if (!is.null(A_eq)) nrow(A_eq) else 0

  b_ub_shifted <- if (rows_ub > 0) b_ub - A_ub %*% shift else numeric(0)
  b_eq_shifted <- if (rows_eq > 0) b_eq - A_eq %*% shift else numeric(0)

  # Upper-bound constraints: x'_i <= ub_i - lb_i  (finite ub only)
  finite_ub <- which(is.finite(ub))
  n_ub_extra <- length(finite_ub)

  # Build combined inequality system (including upper-bound rows)
  if (rows_ub > 0 && n_ub_extra > 0) {
    A_ub_full <- rbind(A_ub,
                       matrix(0, nrow = n_ub_extra, ncol = n))
    for (k in seq_along(finite_ub)) {
      A_ub_full[rows_ub + k, finite_ub[k]] <- 1
    }
    b_ub_full <- c(b_ub_shifted, ub[finite_ub] - lb[finite_ub])
  } else if (rows_ub > 0) {
    A_ub_full <- A_ub
    b_ub_full <- b_ub_shifted
  } else if (n_ub_extra > 0) {
    A_ub_full <- matrix(0, nrow = n_ub_extra, ncol = n)
    for (k in seq_along(finite_ub)) {
      A_ub_full[k, finite_ub[k]] <- 1
    }
    b_ub_full <- ub[finite_ub] - lb[finite_ub]
  } else {
    A_ub_full <- matrix(0, nrow = 0, ncol = n)
    b_ub_full <- numeric(0)
  }

  m_ub <- nrow(A_ub_full)

  # Handle negative RHS in inequalities by multiplying row by -1 (flip to >=)
  # We will treat flipped rows as >= constraints needing surplus + artificial
  flip_idx <- which(b_ub_full < 0)
  ge_from_flip_A <- matrix(0, nrow = 0, ncol = n)
  ge_from_flip_b <- numeric(0)
  if (length(flip_idx) > 0) {
    ge_from_flip_A <- -A_ub_full[flip_idx, , drop = FALSE]
    ge_from_flip_b <- -b_ub_full[flip_idx]
    A_ub_full <- A_ub_full[-flip_idx, , drop = FALSE]
    b_ub_full <- b_ub_full[-flip_idx]
  }
  m_ub_pos <- nrow(A_ub_full)

  # Equality rows also need artificials (no natural slack)
  m_eq <- if (!is.null(A_eq)) nrow(A_eq) else 0
  m_ge <- nrow(ge_from_flip_A)

  n_slack      <- m_ub_pos          # one slack per <= row
  n_surplus    <- m_ge              # one surplus per >= row
  n_art_ge     <- m_ge              # one artificial per >= row
  n_art_eq     <- m_eq              # one artificial per equality row
  n_artificial <- n_art_ge + n_art_eq

  total_vars <- n + n_slack + n_surplus + n_artificial

  # Variable names
  orig_names   <- paste0("x", seq_len(n))
  slack_names  <- if (n_slack   > 0) paste0("s",  seq_len(n_slack))   else character(0)
  surp_names   <- if (n_surplus > 0) paste0("sp", seq_len(n_surplus)) else character(0)
  art_names    <- if (n_artificial > 0) paste0("a", seq_len(n_artificial)) else character(0)
  var_names    <- c(orig_names, slack_names, surp_names, art_names)

  # Total rows
  total_rows <- m_ub_pos + m_ge + m_eq

  if (total_rows == 0) {
    stop("No constraints provided.")
  }

  A_std <- matrix(0, nrow = total_rows, ncol = total_vars)
  b_std <- numeric(total_rows)

  row_ptr <- 1

  # <= rows: add slack
  if (m_ub_pos > 0) {
    rows  <- row_ptr:(row_ptr + m_ub_pos - 1)
    A_std[rows, 1:n] <- A_ub_full
    for (k in seq_len(m_ub_pos)) {
      A_std[rows[k], n + k] <- 1
    }
    b_std[rows] <- b_ub_full
    row_ptr <- row_ptr + m_ub_pos
  }

  # >= rows (from flipped negatives): subtract surplus, add artificial
  if (m_ge > 0) {
    rows <- row_ptr:(row_ptr + m_ge - 1)
    A_std[rows, 1:n] <- ge_from_flip_A
    for (k in seq_len(m_ge)) {
      A_std[rows[k], n + n_slack + k]              <- -1  # surplus
      A_std[rows[k], n + n_slack + n_surplus + k]  <-  1  # artificial
    }
    b_std[rows] <- ge_from_flip_b
    row_ptr <- row_ptr + m_ge
  }

  # Equality rows: add artificial only
  if (m_eq > 0) {
    rows <- row_ptr:(row_ptr + m_eq - 1)
    A_std[rows, 1:n] <- A_eq
    for (k in seq_len(m_eq)) {
      A_std[rows[k], n + n_slack + n_surplus + n_art_ge + k] <- 1
    }
    b_std[rows] <- b_eq_shifted
    row_ptr <- row_ptr + m_eq
  }

  # Objective: artificials get 0 in standard-form objective
  c_std <- c(c_vec, rep(0, n_slack + n_surplus + n_artificial))

  list(
    c              = c_std,
    A              = A_std,
    b              = b_std,
    n_original     = n,
    n_slack        = n_slack,
    n_surplus      = n_surplus,
    n_artificial   = n_artificial,
    variable_names = var_names,
    maximize       = maximize
  )
}
