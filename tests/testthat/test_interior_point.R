library(testthat)
library(lpSolveR)

# ── interior_point ────────────────────────────────────────────────────────────

test_that("returns optimal status on a simple LP", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_equal(res$status, "optimal")
})

test_that("objective value is a finite number", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_true(is.finite(res$obj_value))
})

test_that("matches simplex_method within tolerance on lp_example", {
  data(lp_example, package = "lpSolveR")
  res_ip <- interior_point(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                           maximize = lp_example$maximize)
  res_sp <- simplex_method(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                           maximize = lp_example$maximize)
  
  expect_equal(res_ip$status, "optimal")
  expect_lt(abs(res_ip$obj_value - res_sp$obj_value), 1e-4)
})

test_that("duality_gap_history is a non-empty numeric vector", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_true(is.numeric(res$duality_gap_history))
  expect_true(length(res$duality_gap_history) > 0)
})

test_that("duality gap decreases (is converging)", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  gaps <- res$duality_gap_history
  # The final gap should be less than the initial gap
  expect_lt(tail(gaps, 1), gaps[1])
})

test_that("final duality gap is very small when optimal", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_lt(tail(res$duality_gap_history, 1), 1e-6)
})

test_that("solution vector has the correct length", {
  obj  <- c(1, 2, 3)
  A_ub <- matrix(c(1, 1, 1, 2, 1, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(10, 14)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_equal(length(res$solution), 3)
})

test_that("solution values are non-negative", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_true(all(res$solution >= -1e-6))
})

test_that("returns all required fields", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_true(!is.null(res$solution))
  expect_true(!is.null(res$obj_value))
  expect_true(!is.null(res$status))
  expect_true(!is.null(res$iterations))
  expect_true(!is.null(res$runtime))
  expect_true(!is.null(res$duality_gap_history))
})

test_that("iterations count matches length of duality_gap_history", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_equal(res$iterations, length(res$duality_gap_history))
})

test_that("runtime is non-negative", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)
  
  expect_true(res$runtime >= 0)
})

test_that("works with equality constraints", {
  res <- interior_point(c(1, 2),
                        A_eq = matrix(c(1, 1), nrow = 1),
                        b_eq = 3)
  
  expect_true(res$status %in% c("optimal", "max_iter_reached"))
  expect_true(is.numeric(res$obj_value))
})