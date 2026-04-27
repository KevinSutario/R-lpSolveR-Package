library(testthat)
library(lpSolveR)

# ------------------------------------------------------------------ 1. Standard form
test_that("lp_to_standard_form produces correct dimensions", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)

  sf <- lp_to_standard_form(obj, A_ub, b_ub)

  # 2 original vars + 2 slacks = 4 total columns
  expect_equal(sf$n_original,   2)
  expect_equal(sf$n_slack,      2)
  expect_equal(sf$n_artificial, 0)
  expect_equal(ncol(sf$A), 4)
  expect_equal(nrow(sf$A), 2)
  expect_equal(length(sf$b), 2)
  expect_equal(length(sf$c), 4)
})

test_that("lp_to_standard_form handles equality constraints", {
  obj  <- c(1, 2)
  A_eq <- matrix(c(1, 1), nrow = 1)
  b_eq <- c(3)

  sf <- lp_to_standard_form(obj, A_eq = A_eq, b_eq = b_eq)

  expect_equal(sf$n_slack,      0)
  expect_equal(sf$n_artificial, 1)
  expect_equal(ncol(sf$A), 3)  # 2 orig + 1 artificial
  expect_equal(nrow(sf$A), 1)
})

test_that("lp_to_standard_form flips objective for maximisation", {
  obj <- c(3, 5)
  sf  <- lp_to_standard_form(obj, maximize = TRUE,
                               A_ub = matrix(c(1,0,0,1), nrow=2, byrow=TRUE),
                               b_ub = c(4, 6))
  expect_equal(sf$c[1:2], -obj)
})

# ------------------------------------------------------------------ 2. Simplex
test_that("simplex_method solves a simple 2-variable LP correctly", {
  # max 2x1 + 3x2  s.t.  x1+x2<=4, 2x1+x2<=6, x2<=2,  x1,x2>=0
  # Optimal: x1=2, x2=2 => obj=10
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1,
                   2, 1,
                   0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)

  res <- simplex_method(obj, A_ub, b_ub, maximize = TRUE)

  expect_equal(res$status, "optimal")
  expect_lt(abs(res$obj_value - 10), 1e-6)
  expect_length(res$solution, 2)
})

test_that("simplex_method returns correct fields", {
  obj  <- c(1, 1)
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(3, 3)
  res  <- simplex_method(obj, A_ub, b_ub)

  expect_true(!is.null(res$iterations))
  expect_true(!is.null(res$runtime))
  expect_true(!is.null(res$basis))
  expect_true(!is.null(res$tableau_history))
  expect_true(is.character(res$status))
})

# ------------------------------------------------------------------ 3. Two-phase simplex
test_that("two_phase_simplex detects infeasibility", {
  # x1 >= 1 AND x1 <= 0  -> infeasible
  obj  <- c(1, 0)
  A_ub <- matrix(c(1, 0, -1, 0), nrow = 2, byrow = TRUE)
  b_ub <- c(-1, -1)   # forces x1 >= 1 AND x1 >= 1 after flip, but...
  # Simpler infeasible: Ax=b with b<0 after standard form forces infeasibility
  # Use: x1 + x2 = 5 AND x1 + x2 = 3
  A_eq2 <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  b_eq2 <- c(5, 3)

  res <- two_phase_simplex(c(1, 1), A_eq = A_eq2, b_eq = b_eq2)
  expect_equal(res$status, "infeasible")
})

test_that("two_phase_simplex solves an equality-constrained LP", {
  # min x1 + 2*x2  s.t.  x1 + x2 = 3, x1,x2 >= 0
  # Optimal at x1=3, x2=0 => obj=3
  obj  <- c(1, 2)
  A_eq <- matrix(c(1, 1), nrow = 1)
  b_eq <- c(3)

  res <- two_phase_simplex(obj, A_eq = A_eq, b_eq = b_eq)

  expect_equal(res$status, "optimal")
  expect_lt(abs(res$obj_value - 3), 1e-6)
  expect_true(!is.null(res$phase1_iterations))
  expect_true(!is.null(res$phase2_iterations))
})

# ------------------------------------------------------------------ 4. Interior point
test_that("interior_point solution matches simplex within tolerance", {
  data(lp_example, package = "lpSolveR")

  res_ip  <- interior_point(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                              maximize = lp_example$maximize)
  res_sp  <- simplex_method(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                             maximize = lp_example$maximize)

  expect_equal(res_ip$status, "optimal")
  expect_lt(abs(res_ip$obj_value - res_sp$obj_value), 1e-4)
})

test_that("interior_point returns duality_gap_history", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  res  <- interior_point(obj, A_ub, b_ub)

  expect_true(length(res$duality_gap_history) > 0)
  expect_true(is.numeric(res$duality_gap_history))
})

# ------------------------------------------------------------------ 5. compare_lp_methods
test_that("compare_lp_methods returns summary_df with 3 rows and correct columns", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)

  cmp <- compare_lp_methods(obj, A_ub, b_ub)

  expect_true(is.list(cmp))
  expect_true("summary_df" %in% names(cmp))
  expect_true("results" %in% names(cmp))

  df <- cmp$summary_df
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_true(all(c("method", "obj_value", "iterations",
                     "runtime_sec", "status") %in% names(df)))
})

# ------------------------------------------------------------------ 6. check_lp_feasibility
test_that("check_lp_feasibility returns TRUE for the example dataset", {
  data(lp_example, package = "lpSolveR")

  result <- check_lp_feasibility(lp_example$A_ub, lp_example$b_ub)

  expect_true(is.list(result))
  expect_true(result$feasible)
  expect_true(is.character(result$message))
})

test_that("check_lp_feasibility returns FALSE for contradictory constraints", {
  A_eq <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  b_eq <- c(5, 3)

  result <- check_lp_feasibility(A_eq = A_eq, b_eq = b_eq)
  expect_false(result$feasible)
})
