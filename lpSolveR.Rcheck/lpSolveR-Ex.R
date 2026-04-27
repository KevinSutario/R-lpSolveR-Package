pkgname <- "lpSolveR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('lpSolveR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("check_lp_feasibility")
### * check_lp_feasibility

flush(stderr()); flush(stdout())

### Name: check_lp_feasibility
### Title: Quick Feasibility Check for an LP
### Aliases: check_lp_feasibility

### ** Examples

A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
b_ub <- c(14, 14, 8)
check_lp_feasibility(A_ub, b_ub)




cleanEx()
nameEx("compare_lp_methods")
### * compare_lp_methods

flush(stderr()); flush(stdout())

### Name: compare_lp_methods
### Title: Compare LP Solvers on the Same Problem
### Aliases: compare_lp_methods

### ** Examples

obj  <- c(2, 3)
A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub <- c(4, 6)
cmp  <- compare_lp_methods(obj, A_ub, b_ub)
cmp$summary_df




cleanEx()
nameEx("interior_point")
### * interior_point

flush(stderr()); flush(stdout())

### Name: interior_point
### Title: Solve an LP Using the Primal-Dual Interior Point Method
### Aliases: interior_point

### ** Examples

obj  <- c(2, 3)
A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub <- c(4, 6)
res  <- interior_point(obj, A_ub, b_ub)
res$obj_value




cleanEx()
nameEx("lp_example")
### * lp_example

flush(stderr()); flush(stdout())

### Name: lp_example
### Title: Example LP Dataset
### Aliases: lp_example
### Keywords: datasets

### ** Examples

data(lp_example)
lp_example$obj
lp_example$A_ub




cleanEx()
nameEx("lp_to_standard_form")
### * lp_to_standard_form

flush(stderr()); flush(stdout())

### Name: lp_to_standard_form
### Title: Convert a General LP to Standard Form
### Aliases: lp_to_standard_form

### ** Examples

obj   <- c(2, 3)
A_ub  <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub  <- c(4, 6)
sf    <- lp_to_standard_form(obj, A_ub, b_ub)
sf$A
sf$variable_names




cleanEx()
nameEx("plot_convergence")
### * plot_convergence

flush(stderr()); flush(stdout())

### Name: plot_convergence
### Title: Plot Convergence of LP Solvers
### Aliases: plot_convergence

### ** Examples

obj  <- c(2, 3)
A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub <- c(4, 6)
cmp  <- compare_lp_methods(obj, A_ub, b_ub)
p    <- plot_convergence(cmp)
print(p)




cleanEx()
nameEx("plot_lp_feasible_region")
### * plot_lp_feasible_region

flush(stderr()); flush(stdout())

### Name: plot_lp_feasible_region
### Title: Plot the Feasible Region of a 2-Variable LP
### Aliases: plot_lp_feasible_region

### ** Examples

obj  <- c(2, 3)
A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub <- c(4, 6)
p    <- plot_lp_feasible_region(obj, A_ub, b_ub)
print(p)




cleanEx()
nameEx("simplex_method")
### * simplex_method

flush(stderr()); flush(stdout())

### Name: simplex_method
### Title: Solve an LP Using the Simplex Method
### Aliases: simplex_method

### ** Examples

obj  <- c(2, 3)
A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
b_ub <- c(4, 6)
res  <- simplex_method(obj, A_ub, b_ub)
res$obj_value
res$solution




cleanEx()
nameEx("summarize_lp")
### * summarize_lp

flush(stderr()); flush(stdout())

### Name: summarize_lp
### Title: Print a Human-Readable Summary of an LP Problem
### Aliases: summarize_lp

### ** Examples

obj  <- c(2, 3, 1)
A_ub <- matrix(c(1, 2, 1, 2, 1, 3, 1, 1, 1), nrow = 3, byrow = TRUE)
b_ub <- c(14, 14, 8)
summarize_lp(obj, A_ub, b_ub)




cleanEx()
nameEx("two_phase_simplex")
### * two_phase_simplex

flush(stderr()); flush(stdout())

### Name: two_phase_simplex
### Title: Solve an LP Using the Two-Phase Simplex Method
### Aliases: two_phase_simplex

### ** Examples

obj  <- c(1, 2)
A_eq <- matrix(c(1, 1), nrow = 1)
b_eq <- 3
res  <- two_phase_simplex(obj, A_eq = A_eq, b_eq = b_eq)
res$obj_value




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
