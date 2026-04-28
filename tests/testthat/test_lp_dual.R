library(testthat)
library(lpSolveR)

# ── lp_dual ───────────────────────────────────────────────────────────────────

test_that("constructs dual with correct number of variables", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  dual <- lp_dual(obj, A_ub, b_ub)
  
  # Dual should have one variable per primal constraint
  expect_equal(dual$n_dual_vars, 2)
  expect_equal(length(dual$dual_obj), 2)
})

test_that("dual objective has correct signs for minimization", {
  obj  <- c(2, 3, 1)
  A_ub <- matrix(c(1, 2, 1, 2, 1, 3), nrow = 2, byrow = TRUE)
  b_ub <- c(10, 15)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  # For minimization, dual objective should equal primal RHS (positive)
  expect_equal(dual$dual_obj, b_ub)
})

test_that("dual objective has correct signs for maximization", {
  obj  <- c(2, 3, 1)
  A_ub <- matrix(c(1, 2, 1, 2, 1, 3), nrow = 2, byrow = TRUE)
  b_ub <- c(10, 15)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = TRUE)
  
  # For maximization, dual objective should be negated
  expect_equal(dual$dual_obj, -b_ub)
})

test_that("dual constraint matrix is transpose of primal for minimization", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  # For min problem: dual has A^T y <= c
  expect_equal(dual$dual_A_ub, t(A_ub))
})

test_that("dual RHS equals primal objective for minimization", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  expect_equal(dual$dual_b_ub, obj)
})

test_that("primal min becomes dual max", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  expect_true(dual$dual_maximize)
})

test_that("primal max becomes dual min", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_false(dual$dual_maximize)
})

test_that("handles mixed inequality and equality constraints", {
  obj  <- c(1, 2, 3)
  A_ub <- matrix(c(1, 0, 0), nrow = 1)
  b_ub <- c(5)
  A_eq <- matrix(c(0, 1, 1), nrow = 1)
  b_eq <- c(3)
  
  dual <- lp_dual(obj, A_ub, b_ub, A_eq, b_eq)
  
  # Should have 2 dual variables (1 from inequality, 1 from equality)
  expect_equal(dual$n_dual_vars, 2)
})

test_that("dual construction is internally consistent", {
  # Test that the dual structure is correct
  obj  <- c(3, 5)
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = TRUE)
  
  # Check dimensions
  expect_equal(nrow(dual$dual_A_ub), 2)  # 2 primal variables
  expect_equal(ncol(dual$dual_A_ub), 2)  # 2 dual variables
  expect_equal(length(dual$dual_b_ub), 2)
  expect_equal(length(dual$dual_obj), 2)
})

test_that("returns all required fields", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  
  dual <- lp_dual(obj, A_ub, b_ub)
  
  expect_true(!is.null(dual$dual_obj))
  expect_true(!is.null(dual$dual_A_ub))
  expect_true(!is.null(dual$dual_b_ub))
  expect_true(!is.null(dual$dual_maximize))
  expect_true(!is.null(dual$primal_form))
  expect_true(!is.null(dual$n_dual_vars))
})

test_that("primal_form field is correct", {
  dual_min <- lp_dual(c(1, 2),
                      A_ub = matrix(c(1, 1), nrow = 1),
                      b_ub = 5,
                      maximize = FALSE)
  dual_max <- lp_dual(c(1, 2),
                      A_ub = matrix(c(1, 1), nrow = 1),
                      b_ub = 5,
                      maximize = TRUE)
  
  expect_equal(dual_min$primal_form, "minimisation")
  expect_equal(dual_max$primal_form, "maximisation")
})

test_that("raises error when no constraints provided", {
  expect_error(lp_dual(c(1, 2)))
})

test_that("dual of dual recovers original problem structure", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  # First dual
  dual1 <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  # Dual of dual
  dual2 <- lp_dual(dual1$dual_obj, dual1$dual_A_ub, dual1$dual_b_ub,
                   maximize = dual1$dual_maximize)
  
  # Should have same number of variables as original primal
  expect_equal(dual2$n_dual_vars, length(obj))
})

test_that("dual construction follows theoretical properties", {
  # A well-formed example where we know the relationship
  # Primal: min c^T x s.t. Ax <= b, x >= 0
  # Dual: max b^T y s.t. A^T y <= c, y >= 0
  
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = FALSE)
  
  # Dual should maximize
  expect_true(dual$dual_maximize)
  
  # Dual objective should be primal RHS
  expect_equal(as.vector(dual$dual_obj), b_ub)
  
  # Dual constraints: A^T y <= c means transpose of A
  expect_equal(dual$dual_A_ub, t(A_ub))
  
  # Dual RHS should be primal objective
  expect_equal(as.vector(dual$dual_b_ub), obj)
})

test_that("maximization dual construction follows theoretical properties", {
  # Primal: max c^T x s.t. Ax <= b, x >= 0
  # Dual: min b^T y s.t. A^T y >= c, y >= 0
  # Which becomes: min -b^T y s.t. -A^T y <= -c, y >= 0
  
  obj  <- c(3, 5)
  A_ub <- matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(8, 10)
  
  dual <- lp_dual(obj, A_ub, b_ub, maximize = TRUE)
  
  # Dual should minimize
  expect_false(dual$dual_maximize)
  
  # Dual objective should be negative of primal RHS
  expect_equal(as.vector(dual$dual_obj), -b_ub)
  
  # Dual constraints: -A^T y <= -c
  expect_equal(dual$dual_A_ub, -t(A_ub))
  
  # Dual RHS should be negative of primal objective
  expect_equal(as.vector(dual$dual_b_ub), -obj)
})