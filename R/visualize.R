# Suppress R CMD check note about .data pronoun
utils::globalVariables(".data")

#' Plot the Feasible Region of a 2-Variable LP
#'
#' For linear programs with exactly two decision variables, plots the feasible
#' region as a shaded polygon, draws each constraint boundary, overlays
#' iso-objective contour lines, and marks the optimal vertex.
#'
#' @param obj Numeric vector of length 2: objective coefficients.
#' @param A_ub Numeric matrix with 2 columns: <= inequality constraints.
#' @param b_ub Numeric vector: RHS of the <= constraints.
#' @param maximize Logical. Maximise the objective? Default \code{FALSE}.
#' @param xlim Numeric vector of length 2: x-axis limits. Default \code{c(0, 10)}.
#' @param ylim Numeric vector of length 2: y-axis limits. Default \code{c(0, 10)}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' obj  <- c(2, 3)
#' A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' p    <- plot_lp_feasible_region(obj, A_ub, b_ub)
#' print(p)
#'
#' @export
plot_lp_feasible_region <- function(obj, A_ub, b_ub,
                                     maximize = FALSE,
                                     xlim = c(0, 10),
                                     ylim = c(0, 10)) {
  if (length(obj) != 2) stop("plot_lp_feasible_region only supports 2-variable LPs.")
  if (ncol(A_ub) != 2)  stop("A_ub must have exactly 2 columns.")

  # Solve to get optimal point
  res <- simplex_method(obj, A_ub, b_ub, maximize = maximize)
  opt_x <- if (res$status == "optimal") res$solution[1] else NA
  opt_y <- if (res$status == "optimal") res$solution[2] else NA

  # Compute feasible region vertices via constraint intersections
  m <- nrow(A_ub)
  vertices <- .feasible_vertices(A_ub, b_ub, xlim, ylim)

  # Grid for contour lines
  x_seq <- seq(xlim[1], xlim[2], length.out = 200)
  y_seq <- seq(ylim[1], ylim[2], length.out = 200)
  grid  <- expand.grid(x = x_seq, y = y_seq)
  grid$z <- obj[1] * grid$x + obj[2] * grid$y

  # Contour levels
  z_range <- range(grid$z)
  levels  <- seq(z_range[1], z_range[2], length.out = 6)

  p <- ggplot2::ggplot()

  # Feasible region polygon
  if (nrow(vertices) >= 3) {
    # Order vertices for polygon
    cx <- mean(vertices$x)
    cy <- mean(vertices$y)
    vertices$angle <- atan2(vertices$y - cy, vertices$x - cx)
    vertices <- vertices[order(vertices$angle), ]
    p <- p + ggplot2::geom_polygon(
      data = vertices,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "steelblue", alpha = 0.3, colour = NA
    )
  }

  # Constraint lines
  x_line <- seq(xlim[1], xlim[2], length.out = 300)
  for (i in seq_len(m)) {
    a1 <- A_ub[i, 1]; a2 <- A_ub[i, 2]; rhs <- b_ub[i]
    if (abs(a2) > 1e-10) {
      y_line <- (rhs - a1 * x_line) / a2
      line_df <- data.frame(x = x_line, y = y_line,
                            constraint = paste0("C", i))
      p <- p + ggplot2::geom_line(
        data = line_df,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$constraint),
        linewidth = 0.8
      )
    } else if (abs(a1) > 1e-10) {
      xv <- rhs / a1
      p <- p + ggplot2::geom_vline(
        xintercept = xv,
        colour = grDevices::hcl(h = (i - 1) * 30, l = 50, c = 80),
        linewidth = 0.8
      )
    }
  }

  # Contour lines for objective
  p <- p + ggplot2::geom_contour(
    data = grid,
    ggplot2::aes(x = .data$x, y = .data$y, z = .data$z),
    breaks = levels, colour = "grey40", linetype = "dashed", linewidth = 0.4
  )

  # Optimal point
  if (!is.na(opt_x)) {
    opt_df <- data.frame(x = opt_x, y = opt_y)
    p <- p + ggplot2::geom_point(
      data = opt_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      colour = "red", size = 4, shape = 18
    ) + ggplot2::annotate(
      "text", x = opt_x, y = opt_y,
      label = sprintf("Optimal\n(%.2f, %.2f)", opt_x, opt_y),
      vjust = -0.7, colour = "red", size = 3.2
    )
  }

  direction <- if (maximize) "max" else "min"
  obj_str   <- paste0(direction, " ", obj[1], "x\u2081 + ", obj[2], "x\u2082")

  p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(
      title    = "LP Feasible Region",
      subtitle = obj_str,
      x        = "x\u2081", y = "x\u2082",
      colour   = "Constraint"
    ) +
    ggplot2::theme_bw()
}

#' Plot Convergence of LP Solvers
#'
#' Takes output from \code{\link{compare_lp_methods}} and visualises the
#' convergence behaviour of each solver: objective value per iteration for the
#' simplex methods, and duality gap per iteration for the interior point method.
#'
#' @param compare_result A list returned by \code{\link{compare_lp_methods}}.
#'
#' @return A \code{ggplot} object with facets for each solver.
#'
#' @examples
#' obj  <- c(2, 3)
#' A_ub <- matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)
#' b_ub <- c(4, 6)
#' cmp  <- compare_lp_methods(obj, A_ub, b_ub)
#' p    <- plot_convergence(cmp)
#' print(p)
#'
#' @export
plot_convergence <- function(compare_result) {
  results <- compare_result$results

  plot_data <- list()

  # Simplex methods: tableau_history contains obj values per iteration
  for (method_name in c("simplex", "two_phase")) {
    res <- results[[method_name]]
    th  <- res$tableau_history
    if (length(th) > 0) {
      plot_data[[method_name]] <- data.frame(
        method    = method_name,
        iteration = seq_along(th),
        value     = th,
        metric    = "Objective Value"
      )
    }
  }

  # Interior point: duality gap history
  ip_res <- results$interior_point
  dg     <- ip_res$duality_gap_history
  if (length(dg) > 0) {
    plot_data[["interior_point"]] <- data.frame(
      method    = "interior_point",
      iteration = seq_along(dg),
      value     = dg,
      metric    = "Duality Gap"
    )
  }

  all_data <- do.call(rbind, plot_data)

  if (is.null(all_data) || nrow(all_data) == 0) {
    message("No convergence history to plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  ggplot2::ggplot(all_data,
                  ggplot2::aes(x = .data$iteration,
                               y = .data$value,
                               colour = .data$method)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(
      title  = "LP Solver Convergence",
      x      = "Iteration",
      y      = "Value",
      colour = "Method"
    ) +
    ggplot2::theme_bw()
}

# ---- internal ----------------------------------------------------------------

.feasible_vertices <- function(A_ub, b_ub, xlim, ylim) {
  m <- nrow(A_ub)
  # Augment with non-negativity and box constraints
  A_aug <- rbind(A_ub,
                 c(-1,  0),  # x1 >= 0
                 c( 0, -1),  # x2 >= 0
                 c( 1,  0),  # x1 <= xlim[2]
                 c( 0,  1))  # x2 <= ylim[2]
  b_aug <- c(b_ub, 0, 0, xlim[2], ylim[2])
  n_aug <- nrow(A_aug)

  verts <- data.frame(x = numeric(0), y = numeric(0))
  for (i in seq_len(n_aug - 1)) {
    for (j in (i + 1):n_aug) {
      A2 <- A_aug[c(i, j), ]
      b2 <- b_aug[c(i, j)]
      det_val <- A2[1,1]*A2[2,2] - A2[1,2]*A2[2,1]
      if (abs(det_val) < 1e-10) next
      pt <- solve(A2, b2)
      if (pt[1] >= -1e-8 && pt[2] >= -1e-8 &&
          pt[1] <= xlim[2] + 1e-8 && pt[2] <= ylim[2] + 1e-8) {
        # Check feasibility
        lhs <- as.vector(A_ub %*% pt)
        if (all(lhs <= b_ub + 1e-8)) {
          verts <- rbind(verts, data.frame(x = pt[1], y = pt[2]))
        }
      }
    }
  }
  unique(round(verts, 8))
}
