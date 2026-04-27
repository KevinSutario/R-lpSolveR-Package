library(testthat)
library(lpSolveR)

# ── simplex_method ────────────────────────────────────────────────────────────

test_that("solves a basic minimisation LP at the origin", {
  # min 2x1 + 3x2  s.t. x1<=3, x2<=3 => optimal at (0,0) = 0
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(3, 3)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_equal(res$status, "optimal")
  expect_lt(abs(res$obj_value - 0), 1e-6)
})

test_that("solves a 2-variable maximisation LP correctly", {
  # max 2x1 + 3x2  s.t. x1+x2<=4, 2x1+x2<=6, x2<=2
  # Optimal: x1=2, x2=2, obj=10
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  res  <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_equal(res$status, "optimal")
  expect_lt(abs(res$obj_value - 10), 1e-6)
  expect_equal(length(res$solution), 2)
})

test_that("solution vector has the right length", {
  obj  <- c(1, 2, 3)
  A_ub <- matrix(c(1, 1, 1, 2, 1, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(10, 14)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_equal(length(res$solution), 3)
})

test_that("solution values are non-negative", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_true(all(res$solution >= -1e-9))
})

test_that("returns all required fields", {
  obj  <- c(1, 1)
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(3, 3)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_true(!is.null(res$solution))
  expect_true(!is.null(res$obj_value))
  expect_true(!is.null(res$status))
  expect_true(!is.null(res$iterations))
  expect_true(!is.null(res$runtime))
  expect_true(!is.null(res$basis))
  expect_true(!is.null(res$tableau_history))
})

test_that("status is a character string", {
  obj  <- c(1, 1)
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(3, 3)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_true(is.character(res$status))
})

test_that("iterations is a non-negative integer", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  res  <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_true(res$iterations >= 0)
  expect_true(is.numeric(res$iterations))
})

test_that("runtime is non-negative", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  res  <- simplex_method(obj, A_ub, b_ub)
  
  expect_true(res$runtime >= 0)
})

test_that("tableau_history is a numeric vector", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  res  <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_true(is.numeric(res$tableau_history))
})

test_that("detects infeasibility for contradictory equality constraints", {
  # x1 + x2 = 5  AND  x1 + x2 = 3 is impossible
  res <- simplex_method(c(1, 1),
                        A_eq = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
                        b_eq = c(5, 3))
  
  expect_equal(res$status, "infeasible")
})

test_that("bland rule option does not break results", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  
  res_bland  <- simplex_method(obj, A_ub, b_ub, maximize = TRUE, bland = TRUE)
  res_nobland <- simplex_method(obj, A_ub, b_ub, maximize = TRUE, bland = FALSE)
  
  expect_equal(res_bland$status,    "optimal")
  expect_equal(res_nobland$status,  "optimal")
  expect_lt(abs(res_bland$obj_value - res_nobland$obj_value), 1e-6)
})

test_that("example dataset gives an optimal solution", {
  data(lp_example, package = "lpSolveR")
  res <- simplex_method(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                        maximize = lp_example$maximize)
  
  expect_equal(res$status, "optimal")
  expect_true(is.numeric(res$obj_value))
  expect_false(is.na(res$obj_value))
})