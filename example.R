## Simulation and maximum likelihood estimation of a Matern field on a metric
## tree, under the L2 Hausdorff distance h_kappa.
##
## Companion code for "Repairing covariance functions on metric trees through
## a Hausdorff-type distance".

library(spatstat)
source("aux_functions.R")

## ---- the tree of Figure 1: six edges, four leaves -------------------------
Vx <- c(0.50, 0.24, 0.76, 0.04, 0.46, 0.54, 0.96)
Vy <- c(0.08, 0.40, 0.40, 0.78, 0.90, 0.90, 0.78)
Ed <- matrix(c(1,2, 1,3, 2,4, 2,5, 3,6, 3,7), ncol = 2, byrow = TRUE)
net <- linnet(ppp(Vx, Vy, window = owin(c(0,1), c(0,1))), edges = Ed)

## ---- points along the edges, spaced 0.01 ----------------------------------
X <- as.lpp(pointsOnLines(as.psp(net), eps = 0.01), L = net)
H <- h_dist(X)                       # L2 Hausdorff distance
cat(sprintf("n = %d   leaves = %d   max h = %.2f\n",
            npoints(X), length(leaf_id(net)), max(H)))

## ---- one realisation ------------------------------------------------------
set.seed(1)
nu <- 3/2; sigma2 <- 1; range <- 0.6; nugget <- 0.1
y <- sim_field(H, sigma2, nugget, range, nu, z = rnorm(npoints(X)))

## ---- maximum likelihood, nu fixed -----------------------------------------
fit <- fit_mle(y, H, nu)
print(round(rbind(estimate = fit$est, se = fit$se), 3))
