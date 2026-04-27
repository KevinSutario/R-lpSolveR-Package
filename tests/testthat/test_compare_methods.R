library(testthat)
library(lpSolveR)

# ── compare_lp_methods ────────────────────────────────────────────────────────

test_that("returns a list with results and summary_df", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(is.list(cmp))
  expect_true("results"    %in% names(cmp))
  expect_true("summary_df" %in% names(cmp))
})

test_that("summary_df has exactly 3 rows (one per method)", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_equal(nrow(cmp$summary_df), 3)
})

test_that("summary_df contains the correct column names", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expected_cols <- c("method", "obj_value", "iterations", "runtime_sec", "status")
  expect_true(all(expected_cols %in% names(cmp$summary_df)))
})

test_that("summary_df is a data.frame", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_s3_class(cmp$summary_df, "data.frame")
})

test_that("results list contains simplex, two_phase, and interior_point", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true("simplex"        %in% names(cmp$results))
  expect_true("two_phase"      %in% names(cmp$results))
  expect_true("interior_point" %in% names(cmp$results))
})

test_that("all three methods report 'optimal' status on a solvable LP", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_true(all(cmp$summary_df$status == "optimal"))
})

test_that("all three methods return comparable objective values", {
  data(lp_example, package = "lpSolveR")
  cmp <- compare_lp_methods(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                            maximize = lp_example$maximize)
  
  vals <- cmp$summary_df$obj_value
  # All objective values should be within 1e-3 of each other
  expect_lt(max(vals) - min(vals), 1e-3)
})

test_that("runtime_sec values are non-negative", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(all(cmp$summary_df$runtime_sec >= 0))
})

test_that("iterations are non-negative integers", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(all(cmp$summary_df$iterations >= 0))
})

test_that("each method result contains a solution field", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(!is.null(cmp$results$simplex$solution))
  expect_true(!is.null(cmp$results$two_phase$solution))
  expect_true(!is.null(cmp$results$interior_point$solution))
})

test_that("interior_point result contains duality_gap_history", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(length(cmp$results$interior_point$duality_gap_history) > 0)
})

test_that("two_phase result contains phase1_iterations", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  
  expect_true(!is.null(cmp$results$two_phase$phase1_iterations))
})