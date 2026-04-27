library(testthat)
library(lpSolveR)

# ── plot_lp_feasible_region ───────────────────────────────────────────────────

test_that("returns a ggplot object", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  p    <- plot_lp_feasible_region(obj, A_ub, b_ub)
  
  expect_s3_class(p, "ggplot")
})

test_that("works with maximisation", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  p    <- plot_lp_feasible_region(obj, A_ub, b_ub, maximize = TRUE)
  
  expect_s3_class(p, "ggplot")
})

test_that("custom xlim and ylim are accepted", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  p    <- plot_lp_feasible_region(obj, A_ub, b_ub, xlim = c(0, 7), ylim = c(0, 7))
  
  expect_s3_class(p, "ggplot")
})

test_that("raises error when obj has more than 2 elements", {
  obj  <- c(1, 2, 3)
  A_ub <- matrix(c(1, 1, 1), nrow = 1)
  b_ub <- c(5)
  
  expect_error(plot_lp_feasible_region(obj, A_ub, b_ub))
})

test_that("raises error when A_ub does not have exactly 2 columns", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 1, 2, 1, 1), nrow = 2, byrow = TRUE)  # 3 cols
  b_ub <- c(4, 6)
  
  expect_error(plot_lp_feasible_region(obj, A_ub, b_ub))
})

test_that("gg object has expected layers (polygon, line, point)", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  p    <- plot_lp_feasible_region(obj, A_ub, b_ub)
  
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(length(layer_classes) > 0)
})

# ── plot_convergence ──────────────────────────────────────────────────────────

test_that("returns a ggplot object", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  p    <- plot_convergence(cmp)
  
  expect_s3_class(p, "ggplot")
})

test_that("works on lp_example dataset comparison", {
  data(lp_example, package = "lpSolveR")
  cmp <- compare_lp_methods(lp_example$obj, lp_example$A_ub, lp_example$b_ub,
                            maximize = lp_example$maximize)
  p   <- plot_convergence(cmp)
  
  expect_s3_class(p, "ggplot")
})

test_that("plot has at least one layer", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1, 0, 1), nrow = 3, byrow = TRUE)
  b_ub <- c(4, 6, 2)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub, maximize = TRUE)
  p    <- plot_convergence(cmp)
  
  expect_true(length(p$layers) > 0)
})

test_that("plot uses facets (one panel per metric)", {
  obj  <- c(2, 3)
  A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
  b_ub <- c(4, 6)
  cmp  <- compare_lp_methods(obj, A_ub, b_ub)
  p    <- plot_convergence(cmp)
  
  expect_false(is.null(p$facet))
})