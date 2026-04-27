#' Compare LP Solvers on the Same Problem
#'
#' Runs \code{\link{simplex_method}}, \code{\link{two_phase_simplex}}, and
#' \code{\link{interior_point}} on an identical LP and returns a side-by-side
#' comparison of their results and performance.
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
#' @param tol Numeric tolerance passed to each solver. Default \code{1e-8}.
#' @param max_iter Maximum iterations passed to each solver. Default \code{1000}.
#'
#' @return A named list:
#' \describe{
#'   \item{results}{Named list with elements \code{simplex},
#'     \code{two_phase}, and \code{interior_point}, each containing the
#'     full output of the respective solver.}
#'   \item{summary_df}{A \code{data.frame} with columns \code{method},
#'     \code{obj_value}, \code{iterations}, \code{runtime_sec}, and
#'     \code{status}.}
#' }
#'
#' @examples
#' obj  <- c(2, 3)
#' A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' cmp  <- compare_lp_methods(obj, A_ub, b_ub)
#' cmp$summary_df
#'
#' @export
compare_lp_methods <- function(obj,
                                A_ub     = NULL, b_ub = NULL,
                                A_eq     = NULL, b_eq = NULL,
                                maximize = FALSE,
                                tol      = 1e-8,
                                max_iter = 1000) {
  run_timed <- function(fn, ...) {
    t <- system.time(res <- fn(...))
    res$runtime <- as.numeric(t["elapsed"])
    res
  }

  res_simplex <- run_timed(simplex_method,
                           obj, A_ub, b_ub, A_eq, b_eq,
                           maximize = maximize, tol = tol,
                           max_iter = max_iter)

  res_two <- run_timed(two_phase_simplex,
                       obj, A_ub, b_ub, A_eq, b_eq,
                       maximize = maximize, tol = tol,
                       max_iter = max_iter)

  res_ip <- run_timed(interior_point,
                      obj, A_ub, b_ub, A_eq, b_eq,
                      maximize = maximize, tol = tol,
                      max_iter = max_iter)

  results <- list(
    simplex        = res_simplex,
    two_phase      = res_two,
    interior_point = res_ip
  )

  summary_df <- data.frame(
    method      = c("simplex", "two_phase", "interior_point"),
    obj_value   = c(res_simplex$obj_value,
                    res_two$obj_value,
                    res_ip$obj_value),
    iterations  = c(res_simplex$iterations,
                    res_two$iterations,
                    res_ip$iterations),
    runtime_sec = c(res_simplex$runtime,
                    res_two$runtime,
                    res_ip$runtime),
    status      = c(res_simplex$status,
                    res_two$status,
                    res_ip$status),
    stringsAsFactors = FALSE
  )

  list(results = results, summary_df = summary_df)
}
