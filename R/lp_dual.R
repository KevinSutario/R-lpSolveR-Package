#' Construct the Dual LP Problem
#'
#' Given a primal linear program in standard inequality form, constructs the
#' corresponding dual problem. The dual provides bounds on the optimal value
#' and can sometimes be easier to solve than the primal.
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
#' @param maximize Logical. Is the primal a maximisation problem? Default \code{FALSE}.
#'
#' @return A named list containing the dual problem:
#' \describe{
#'   \item{dual_obj}{Objective coefficients for the dual problem.}
#'   \item{dual_A_ub}{Inequality constraint matrix for the dual.}
#'   \item{dual_b_ub}{Inequality RHS vector for the dual.}
#'   \item{dual_A_eq}{Equality constraint matrix for the dual (if any).}
#'   \item{dual_b_eq}{Equality RHS vector for the dual (if any).}
#'   \item{dual_maximize}{Logical indicating dual objective direction.}
#'   \item{primal_form}{Character describing the primal form assumed.}
#'   \item{n_dual_vars}{Number of dual variables.}
#' }
#'
#' @details
#' For a primal minimisation problem of the form:
#' \deqn{\min c^T x \text{ subject to } Ax \leq b, x \geq 0}
#' The dual is:
#' \deqn{\max b^T y \text{ subject to } A^T y \leq c, y \geq 0}
#'
#' For a primal maximisation problem:
#' \deqn{\max c^T x \text{ subject to } Ax \leq b, x \geq 0}
#' The dual is:
#' \deqn{\min b^T y \text{ subject to } A^T y \geq c, y \geq 0}
#'
#' @examples
#' # Primal: max 3x1 + 5x2 subject to x1 <= 4, x2 <= 6
#' obj  <- c(3, 5)
#' A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' dual <- lp_dual(obj, A_ub, b_ub, maximize = TRUE)
#' 
#' # Solve both primal and dual
#' primal_res <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)
#' dual_res   <- simplex_method(dual$dual_obj, dual$dual_A_ub, dual$dual_b_ub,
#'                               maximize = dual$dual_maximize)
#'
#' @export
lp_dual <- function(obj,
                    A_ub = NULL, b_ub = NULL,
                    A_eq = NULL, b_eq = NULL,
                    maximize = FALSE) {
  n <- length(obj)
  
  # Determine problem dimensions
  m_ub <- if (!is.null(A_ub)) nrow(A_ub) else 0
  m_eq <- if (!is.null(A_eq)) nrow(A_eq) else 0
  
  if (m_ub == 0 && m_eq == 0) {
    stop("At least one constraint (inequality or equality) must be provided.")
  }
  
  # Number of dual variables = number of primal constraints
  n_dual <- m_ub + m_eq
  
  # Build constraint matrix A^T
  A_transpose <- matrix(0, nrow = n, ncol = n_dual)
  if (m_ub > 0) {
    A_transpose[, 1:m_ub] <- t(A_ub)
  }
  if (m_eq > 0) {
    A_transpose[, (m_ub + 1):n_dual] <- t(A_eq)
  }
  
  # Dual objective coefficients = primal RHS
  dual_obj <- numeric(n_dual)
  if (m_ub > 0) dual_obj[1:m_ub] <- b_ub
  if (m_eq > 0) dual_obj[(m_ub + 1):n_dual] <- b_eq
  
  if (maximize) {
    # Primal: max c^T x s.t. Ax <= b, x >= 0
    # Dual:   min b^T y s.t. A^T y >= c, y >= 0
    # Negate dual objective for use in minimization
    dual_obj <- -dual_obj
    # Convert A^T y >= c to -A^T y <= -c
    dual_A_ub <- -A_transpose
    dual_b_ub <- -obj
    dual_maximize <- FALSE
  } else {
    # Primal: min c^T x s.t. Ax <= b, x >= 0  
    # Dual:   max b^T y s.t. A^T y <= c, y >= 0
    dual_A_ub <- A_transpose
    dual_b_ub <- obj
    dual_maximize <- TRUE
  }
  
  dual_A_eq <- NULL
  dual_b_eq <- NULL
  
  primal_form <- if (maximize) "maximisation" else "minimisation"
  
  list(
    dual_obj      = dual_obj,
    dual_A_ub     = dual_A_ub,
    dual_b_ub     = dual_b_ub,
    dual_A_eq     = dual_A_eq,
    dual_b_eq     = dual_b_eq,
    dual_maximize = dual_maximize,
    primal_form   = primal_form,
    n_dual_vars   = n_dual
  )
}