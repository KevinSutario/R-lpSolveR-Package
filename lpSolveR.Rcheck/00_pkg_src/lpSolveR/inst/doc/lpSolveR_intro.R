knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)
library(lpSolveR)

data(lp_example)
str(lp_example)

summarize_lp(
  obj      = lp_example$obj,
  A_ub     = lp_example$A_ub,
  b_ub     = lp_example$b_ub,
  maximize = lp_example$maximize
)

sf <- lp_to_standard_form(
  obj      = lp_example$obj,
  A_ub     = lp_example$A_ub,
  b_ub     = lp_example$b_ub,
  maximize = lp_example$maximize
)

cat("Original variables :", sf$n_original, "\n")
cat("Slack variables    :", sf$n_slack,    "\n")
cat("Artificial vars    :", sf$n_artificial, "\n")
cat("Standard-form A dimensions:", nrow(sf$A), "x", ncol(sf$A), "\n")
sf$variable_names

res_simplex <- simplex_method(
  obj      = lp_example$obj,
  A_ub     = lp_example$A_ub,
  b_ub     = lp_example$b_ub,
  maximize = lp_example$maximize
)

cat("Status     :", res_simplex$status, "\n")
cat("Objective  :", res_simplex$obj_value, "\n")
cat("Solution   :", res_simplex$solution, "\n")
cat("Iterations :", res_simplex$iterations, "\n")
cat("Runtime (s):", round(res_simplex$runtime, 4), "\n")

# Equality-constrained LP: min x1 + 2*x2  s.t.  x1 + x2 = 3
res_two <- two_phase_simplex(
  obj  = c(1, 2),
  A_eq = matrix(c(1, 1), nrow = 1),
  b_eq = 3
)

cat("Status            :", res_two$status, "\n")
cat("Objective value   :", res_two$obj_value, "\n")
cat("Solution          :", res_two$solution, "\n")
cat("Phase 1 iterations:", res_two$phase1_iterations, "\n")
cat("Phase 2 iterations:", res_two$phase2_iterations, "\n")

# Contradictory equalities: x1+x2=5 AND x1+x2=3
res_infeas <- two_phase_simplex(
  obj  = c(1, 1),
  A_eq = matrix(c(1,1,1,1), nrow = 2, byrow = TRUE),
  b_eq = c(5, 3)
)
cat("Status:", res_infeas$status, "\n")

res_ip <- interior_point(
  obj      = lp_example$obj,
  A_ub     = lp_example$A_ub,
  b_ub     = lp_example$b_ub,
  maximize = lp_example$maximize
)

cat("Status     :", res_ip$status, "\n")
cat("Objective  :", round(res_ip$obj_value, 6), "\n")
cat("Solution   :", round(res_ip$solution, 6), "\n")
cat("Iterations :", res_ip$iterations, "\n")
cat("Final gap  :", round(tail(res_ip$duality_gap_history, 1), 8), "\n")

cmp <- compare_lp_methods(
  obj      = lp_example$obj,
  A_ub     = lp_example$A_ub,
  b_ub     = lp_example$b_ub,
  maximize = lp_example$maximize
)

knitr::kable(cmp$summary_df, digits = 6,
             caption = "Side-by-side comparison of LP solvers")

p_conv <- plot_convergence(cmp)
print(p_conv)

obj2  <- c(2, 3)
A_ub2 <- matrix(c(1, 1,
                   2, 1,
                   0, 1), nrow = 3, byrow = TRUE)
b_ub2 <- c(4, 6, 3)

p_region <- plot_lp_feasible_region(obj2, A_ub2, b_ub2,
                                     xlim = c(0, 7), ylim = c(0, 7))
print(p_region)
