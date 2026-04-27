library(testthat)
library(lpSolveR)

# ── lp_to_standard_form ──────────────────────────────────────────────────────

test_that("basic inequality LP: correct variable counts and dimensions", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  sf   <- lp_to_standard_form(obj, A_ub, b_ub)
  
  expect_equal(sf$n_original,   2)
  expect_equal(sf$n_slack,      2)
  expect_equal(sf$n_artificial, 0)
  expect_equal(nrow(sf$A), 2)
  expect_equal(ncol(sf$A), 4)   # 2 orig + 2 slack
  expect_equal(length(sf$b), 2)
  expect_equal(length(sf$c), 4)
})

test_that("slack columns form an identity sub-matrix", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(3, 4, 5)
  sf   <- lp_to_standard_form(obj, A_ub, b_ub)
  
  # Slack columns (3, 4, 5) should be an identity matrix
  slack_cols <- sf$A[, 3:5]
  expect_equal(slack_cols, diag(3))
})

test_that("equality constraint adds one artificial variable, no slacks", {
  obj  <- c(1, 2)
  A_eq <- matrix(c(1, 1), nrow = 1)
  b_eq <- 3
  sf   <- lp_to_standard_form(obj, A_eq = A_eq, b_eq = b_eq)
  
  expect_equal(sf$n_slack,      0)
  expect_equal(sf$n_artificial, 1)
  expect_equal(ncol(sf$A), 3)   # 2 orig + 1 artificial
  expect_equal(nrow(sf$A), 1)
})

test_that("maximisation flips objective sign", {
  obj <- c(3, 5)
  sf  <- lp_to_standard_form(obj, maximize = TRUE,
                             A_ub = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE),
                             b_ub = c(4, 6))
  expect_equal(sf$c[1:2], -obj)
})

test_that("minimisation leaves objective sign unchanged", {
  obj <- c(3, 5)
  sf  <- lp_to_standard_form(obj, maximize = FALSE,
                             A_ub = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE),
                             b_ub = c(4, 6))
  expect_equal(sf$c[1:2], obj)
})

test_that("variable names are correctly generated", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  sf   <- lp_to_standard_form(obj, A_ub, b_ub)
  
  expect_equal(sf$variable_names, c("x1", "x2", "s1"))
})

test_that("variable names include artificials for equality constraints", {
  obj  <- c(1, 2)
  A_eq <- matrix(c(1, 1), nrow = 1)
  b_eq <- 3
  sf   <- lp_to_standard_form(obj, A_eq = A_eq, b_eq = b_eq)
  
  expect_true("a1" %in% sf$variable_names)
})

test_that("RHS vector b is non-negative after standard form conversion", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  sf   <- lp_to_standard_form(obj, A_ub, b_ub)
  
  expect_true(all(sf$b >= 0))
})

test_that("mixed inequality and equality constraints are handled", {
  obj  <- c(1, 2, 3)
  A_ub <- matrix(c(1, 0, 0), nrow = 1)
  b_ub <- c(5)
  A_eq <- matrix(c(0, 1, 1), nrow = 1)
  b_eq <- c(3)
  sf   <- lp_to_standard_form(obj, A_ub, b_ub, A_eq, b_eq)
  
  expect_equal(sf$n_slack,      1)
  expect_equal(sf$n_artificial, 1)
  expect_equal(nrow(sf$A), 2)
})

test_that("error is raised when no constraints are provided", {
  expect_error(lp_to_standard_form(c(1, 2)))
})

test_that("maximize field is recorded in output", {
  sf_min <- lp_to_standard_form(c(1, 2),
                                A_ub = matrix(c(1, 1), nrow = 1), b_ub = 5,
                                maximize = FALSE)
  sf_max <- lp_to_standard_form(c(1, 2),
                                A_ub = matrix(c(1, 1), nrow = 1), b_ub = 5,
                                maximize = TRUE)
  
  expect_false(sf_min$maximize)
  expect_true(sf_max$maximize)
})