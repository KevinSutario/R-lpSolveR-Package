#' Print a Human-Readable Summary of an LP Problem
#'
#' Displays the number of variables, constraints by type, objective direction,
#' and the LP formulation in a formatted layout. Returns invisibly.
#'
#' @param obj Numeric vector of objective coefficients (length n).
#' @param A_ub Numeric matrix of <= inequality constraint coefficients.
#'   Default \code{NULL}.
#' @param b_ub Numeric vector of <= inequality RHS values.
#'   Default \code{NULL}.
#' @param A_eq Numeric matrix of equality constraint coefficients.
#'   Default \code{NULL}.
#' @param b_eq Numeric vector of equality RHS values.
#'   Default \code{NULL}.
#' @param lb Numeric vector of variable lower bounds. Default \code{NULL}.
#' @param ub Numeric vector of variable upper bounds. Default \code{NULL}.
#' @param maximize Logical. Is the objective to be maximised? Default \code{FALSE}.
#'
#' @return Invisibly returns a character string of the summary.
#'
#' @examples
#' obj  <- c(2, 3, 1)
#' A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
#' b_ub <- c(14, 14, 8)
#' summarize_lp(obj, A_ub, b_ub)
#'
#' @export
summarize_lp <- function(obj,
                          A_ub = NULL, b_ub = NULL,
                          A_eq = NULL, b_eq = NULL,
                          lb   = NULL, ub   = NULL,
                          maximize = FALSE) {
  n      <- length(obj)
  m_ub   <- if (!is.null(A_ub)) nrow(A_ub) else 0
  m_eq   <- if (!is.null(A_eq)) nrow(A_eq) else 0
  dir    <- if (maximize) "max" else "min"

  .fmt_expr <- function(coefs) {
    terms <- character(0)
    for (i in seq_along(coefs)) {
      coef <- coefs[i]
      if (coef == 0) next
      sign_str <- if (length(terms) == 0) {
        if (coef < 0) "-" else ""
      } else {
        if (coef < 0) " - " else " + "
      }
      abs_c <- abs(coef)
      c_str <- if (abs_c == 1) "" else as.character(abs_c)
      terms <- c(terms, paste0(sign_str, c_str, "x", i))
    }
    if (length(terms) == 0) "0" else paste(terms, collapse = "")
  }

  lines <- character(0)
  lines <- c(lines, "=== LP Problem Summary ===")
  lines <- c(lines, sprintf("Variables          : %d", n))
  lines <- c(lines, sprintf("Inequality (<=)    : %d", m_ub))
  lines <- c(lines, sprintf("Equality           : %d", m_eq))
  lines <- c(lines, sprintf("Objective direction: %s", toupper(dir)))
  lines <- c(lines, "")
  lines <- c(lines, sprintf("%s  %s", dir, .fmt_expr(obj)))
  lines <- c(lines, "subject to:")

  if (m_ub > 0) {
    for (i in seq_len(m_ub)) {
      lines <- c(lines, sprintf("  %s  <=  %g", .fmt_expr(A_ub[i, ]), b_ub[i]))
    }
  }
  if (m_eq > 0) {
    for (i in seq_len(m_eq)) {
      lines <- c(lines, sprintf("  %s  =  %g", .fmt_expr(A_eq[i, ]), b_eq[i]))
    }
  }

  # Bounds
  lb_eff <- if (!is.null(lb)) lb else rep(0, n)
  ub_eff <- if (!is.null(ub)) ub else rep(Inf, n)
  bound_parts <- character(n)
  for (i in seq_len(n)) {
    lo <- lb_eff[i]; hi <- ub_eff[i]
    if (lo == 0 && is.infinite(hi)) {
      bound_parts[i] <- sprintf("x%d >= 0", i)
    } else if (is.infinite(hi)) {
      bound_parts[i] <- sprintf("x%d >= %g", i, lo)
    } else {
      bound_parts[i] <- sprintf("%g <= x%d <= %g", lo, i, hi)
    }
  }
  lines <- c(lines, paste0("  ", paste(bound_parts, collapse = ", ")))
  lines <- c(lines, "==========================")

  out <- paste(lines, collapse = "\n")
  cat(out, "\n")
  invisible(out)
}


#' Quick Feasibility Check for an LP
#'
#' Runs Phase 1 of the two-phase simplex method to determine whether a linear
#' program has a feasible solution, without solving the full problem.
#'
#' @param A_ub Numeric matrix of <= inequality constraint coefficients.
#'   Default \code{NULL}.
#' @param b_ub Numeric vector of <= inequality RHS values.
#'   Default \code{NULL}.
#' @param A_eq Numeric matrix of equality constraint coefficients.
#'   Default \code{NULL}.
#' @param b_eq Numeric vector of equality RHS values.
#'   Default \code{NULL}.
#'
#' @return A named list:
#' \describe{
#'   \item{feasible}{Logical. \code{TRUE} if a feasible solution exists.}
#'   \item{message}{A character string describing the result.}
#' }
#'
#' @examples
#' A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
#' b_ub <- c(14, 14, 8)
#' check_lp_feasibility(A_ub, b_ub)
#'
#' @export
check_lp_feasibility <- function(A_ub = NULL, b_ub = NULL,
                                  A_eq = NULL, b_eq = NULL) {
  n <- if (!is.null(A_ub)) ncol(A_ub) else if (!is.null(A_eq)) ncol(A_eq) else
    stop("At least one of A_ub or A_eq must be provided.")

  # Use a zero objective — we only care about Phase 1
  obj_dummy <- rep(0, n)

  res <- two_phase_simplex(obj_dummy, A_ub, b_ub, A_eq, b_eq,
                            maximize = FALSE, tol = 1e-8, max_iter = 500)

  if (res$status == "infeasible") {
    list(feasible = FALSE,
         message  = "Problem is infeasible: no feasible solution exists.")
  } else if (res$status %in% c("optimal", "unbounded")) {
    list(feasible = TRUE,
         message  = "Problem is feasible: a basic feasible solution was found.")
  } else {
    list(feasible = NA,
         message  = paste("Feasibility check inconclusive:", res$status))
  }
}
