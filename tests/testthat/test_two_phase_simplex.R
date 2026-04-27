library(testthat)
library(lpSolveR)

# ── two_phase_simplex ─────────────────────────────────────────────────────────

test_that("solves a simple equality-constrained LP", {
  # min x1 + 2*x2  s.t. x1 + x2 = 3 => optimal at x1=3, x2=0, obj=3
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_equal(res$status, "optimal")
  expect_lt(abs(res$obj_value - 3), 1e-6)
})

test_that("detects infeasibility for contradictory equalities", {
  # x1+x2=5 AND x1+x2=3 has no solution
  res <- two_phase_simplex(c(1, 1),
                           A_eq = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
                           b_eq = c(5, 3))
  
  expect_equal(res$status, "infeasible")
})

test_that("solution vector has the correct length", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_equal(length(res$solution), 2)
})

test_that("solution values are non-negative", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_true(all(res$solution >= -1e-9))
})

test_that("returns phase1_iterations and phase2_iterations", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_true(!is.null(res$phase1_iterations))
  expect_true(!is.null(res$phase2_iterations))
  expect_true(res$phase1_iterations >= 0)
  expect_true(res$phase2_iterations >= 0)
})

test_that("total iterations equals phase1 + phase2", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_equal(res$iterations, res$phase1_iterations + res$phase2_iterations)
})

test_that("returns all required fields", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_true(!is.null(res$solution))
  expect_true(!is.null(res$obj_value))
  expect_true(!is.null(res$status))
  expect_true(!is.null(res$iterations))
  expect_true(!is.null(res$runtime))
  expect_true(!is.null(res$basis))
  expect_true(!is.null(res$tableau_history))
  expect_true(!is.null(res$phase1_iterations))
  expect_true(!is.null(res$phase2_iterations))
})

test_that("handles inequality-only LP (no equality constraints)", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- two_phase_simplex(obj, A_ub, b_ub)
  
  expect_equal(res$status, "optimal")
  expect_true(is.numeric(res$obj_value))
})

test_that("matches simplex_method on a standard inequality LP", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  
  res_sp  <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)
  res_tp  <- two_phase_simplex(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_equal(res_sp$status, "optimal")
  expect_equal(res_tp$status, "optimal")
  expect_lt(abs(res_sp$obj_value - res_tp$obj_value), 1e-6)
})

test_that("runtime is non-negative", {
  res <- two_phase_simplex(c(1, 2),
                           A_eq = matrix(c(1, 1), nrow = 1),
                           b_eq = 3)
  
  expect_true(res$runtime >= 0)
})

test_that("handles the lp_example dataset", {
  data(lp_example, package = "lpSolveR")
  res <- two_phase_simplex(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                           maximize = lp_example$maximize)
  
  expect_equal(res$status, "optimal")
  expect_false(is.na(res$obj_value))
})

test_that("obj_value is NA when infeasible", {
  res <- two_phase_simplex(c(1, 1),
                           A_eq = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
                           b_eq = c(5, 3))
  
  expect_true(is.na(res$obj_value))
})