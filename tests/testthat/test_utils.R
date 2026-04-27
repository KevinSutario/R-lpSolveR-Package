library(testthat)
library(lpSolveR)

# ── summarize_lp ──────────────────────────────────────────────────────────────

test_that("returns invisibly a character string", {
  obj  <- c(2, 3, 1)
  A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(14, 14, 8)
  
  out <- capture.output(result <- summarize_lp(obj, A_ub, b_ub))
  expect_true(is.character(result))
})

test_that("output contains variable and constraint counts", {
  obj  <- c(2, 3, 1)
  A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(14, 14, 8)
  
  out <- capture.output(summarize_lp(obj, A_ub, b_ub))
  combined <- paste(out, collapse = "\n")
  
  expect_true(grepl("3", combined))         # 3 variables
  expect_true(grepl("MIN|min", combined))   # minimisation direction
})

test_that("output changes for maximisation", {
  obj  <- c(1, 2)
  A_ub <- matrix(c(1, 1), nrow = 1)
  b_ub <- c(5)
  
  out_min <- capture.output(summarize_lp(obj, A_ub, b_ub, maximize = FALSE))
  out_max <- capture.output(summarize_lp(obj, A_ub, b_ub, maximize = TRUE))
  
  expect_true(grepl("MIN|min", paste(out_min, collapse = "")))
  expect_true(grepl("MAX|max", paste(out_max, collapse = "")))
})

test_that("works when A_eq is provided", {
  obj  <- c(1, 2)
  A_eq <- matrix(c(1, 1), nrow = 1)
  b_eq <- c(3)
  
  expect_silent(capture.output(summarize_lp(obj, A_eq = A_eq, b_eq = b_eq)))
})

test_that("works on lp_example dataset", {
  data(lp_example, package = "lpSolveR")
  expect_silent(capture.output(
    summarize_lp(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                 maximize = lp_example$maximize)
  ))
})

# ── check_lp_feasibility ─────────────────────────────────────────────────────

test_that("returns a list with 'feasible' and 'message' fields", {
  A_ub   <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub   <- c(4, 6)
  result <- check_lp_feasibility(A_ub, b_ub)
  
  expect_true(is.list(result))
  expect_true("feasible" %in% names(result))
  expect_true("message"  %in% names(result))
})

test_that("returns TRUE for the lp_example dataset", {
  data(lp_example, package = "lpSolveR")
  result <- check_lp_feasibility(lp_example$A_ub, lp_example$b_ub)
  
  expect_true(result$feasible)
})

test_that("returns FALSE for contradictory equality constraints", {
  A_eq <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
  b_eq <- c(5, 3)
  
  result <- check_lp_feasibility(A_eq = A_eq, b_eq = b_eq)
  expect_false(result$feasible)
})

test_that("message is a non-empty character string", {
  A_ub   <- matrix(c(1, 1), nrow = 1)
  b_ub   <- c(5)
  result <- check_lp_feasibility(A_ub, b_ub)
  
  expect_true(is.character(result$message))
  expect_true(nchar(result$message) > 0)
})

test_that("feasible field is logical", {
  A_ub   <- matrix(c(1, 1), nrow = 1)
  b_ub   <- c(5)
  result <- check_lp_feasibility(A_ub, b_ub)
  
  expect_true(is.logical(result$feasible))
})

test_that("returns TRUE for a clearly feasible inequality LP", {
  A_ub <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(10, 10)
  
  result <- check_lp_feasibility(A_ub, b_ub)
  expect_true(result$feasible)
})

test_that("raises an error when neither A_ub nor A_eq is provided", {
  expect_error(check_lp_feasibility())
})

test_that("accepts only A_eq without A_ub", {
  A_eq   <- matrix(c(1, 1), nrow = 1)
  b_eq   <- c(3)
  result <- check_lp_feasibility(A_eq = A_eq, b_eq = b_eq)
  
  expect_true(is.list(result))
  expect_true("feasible" %in% names(result))
})