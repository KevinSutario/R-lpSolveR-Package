#' lpSolveR: Linear Programming Solvers in R
#'
#' Provides implementations of the simplex method, two-phase simplex, and
#' primal-dual interior point method for solving linear programming problems.
#' Includes utilities for standard-form conversion, feasibility checking,
#' method comparison, and ggplot2-based visualisation.
#'
#' @docType package
#' @name lpSolveR
"_PACKAGE"

#' Example LP Dataset
#'
#' A small production linear program used in examples and tests throughout the
#' package.  The objective is to minimise cost subject to three resource
#' constraints over three variables.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{obj}{Numeric vector of objective coefficients \code{c(2, 3, 1)}.}
#'   \item{A_ub}{3 x 3 numeric matrix of <= constraint coefficients.}
#'   \item{b_ub}{Numeric vector of RHS values \code{c(14, 14, 8)}.}
#'   \item{maximize}{Logical \code{FALSE} (minimisation problem).}
#'   \item{description}{Character description of the problem.}
#' }
#'
#' @examples
#' data(lp_example)
#' lp_example$obj
#' lp_example$A_ub
#'
#' @source Constructed for illustrative purposes.
"lp_example"
