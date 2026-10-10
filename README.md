# Covariance functions on metric trees

Companion code for *Repairing covariance functions on metric trees through a
Hausdorff-type distance*.

## Files

- `aux_functions.R` — the $L^2$ Hausdorff distance `h_dist` built from the
  distances to the leaves of a tree, the shortest-path distance `d_dist` for
  comparison, the Matern correlation parameterised by its practical range,
  simulation of a Gaussian field from a distance matrix, and maximum
  likelihood estimation.
- `example.R` — a worked example on the tree of Figure 1 of the paper: the
  metric is evaluated at 278 points along the edges, a Matern field with
  smoothness 3/2 is simulated, and its parameters are estimated by maximum
  likelihood.

## Requirements

R with the `spatstat` package, used only to represent the tree as a linear
network and to compute distances along it.
