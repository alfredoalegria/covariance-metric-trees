## leaves of a network: vertices belonging to a single edge
leaf_id <- function(L) {
  deg <- table(c(L$from, L$to))
  as.integer(names(deg)[deg == 1])
}

## shortest-path distance matrix for an lpp
d_dist <- function(X) pairdist(X)

## h distance matrix for an lpp on a tree
h_dist <- function(X, w = NULL) {
  L  <- domain(X)
  lv <- leaf_id(L)
  V  <- vertices(L)
  n  <- npoints(X)
  k  <- length(lv)
  XL <- as.lpp(x = c(coords(X)$x, V$x[lv]),
               y = c(coords(X)$y, V$y[lv]), L = L)
  A  <- pairdist(XL)[1:n, n + seq_len(k), drop = FALSE]
  if (is.null(w)) w <- rep(1/k, k)
  G  <- tcrossprod(A * rep(sqrt(w), each = n))
  dg <- diag(G)
  sqrt(pmax(outer(dg, dg, "+") - 2 * G, 0))
}

#---------------------------------------------------------------------------
## constant c(nu) such that the Matern correlation equals 0.05 at r = range
mat_const <- function(nu)
  uniroot(function(x) 2^(1-nu)/gamma(nu) * x^nu * besselK(x, nu) - 0.05,
          c(1e-8, 1e4))$root

## Matern correlation
mat_cor <- function(r, range, nu) {
  x <- mat_const(nu) * r / range
  out <- 2^(1-nu)/gamma(nu) * x^nu * besselK(x, nu)
  out[r == 0] <- 1
  out
}

## covariance matrix from a distance matrix
cov_mat <- function(D, sigma2, nugget, range, nu)
  sigma2 * mat_cor(D, range, nu) + nugget * diag(nrow(D))

#---------------------------------------------------------------------------

## one realisation from a given covariance matrix, driven by a supplied noise
sim_field <- function(D, sigma2, nugget, range, nu, z, mu = 0) {
  K <- cov_mat(D, sigma2, nugget, range, nu) + diag(1e-9, nrow(D))
  as.vector(mu + t(chol(K)) %*% z)
}

#---------------------------------------------------------------------------

## par = (mu, log sigma2, log nugget, log range); nu fixed
loglik <- function(par, datos, D, nu) {
  lik <- 1.0e8

  mm <- cov_mat(D, exp(par[2]), exp(par[3]), exp(par[4]), nu)
  cholm <- try(chol(mm), silent = TRUE)
  if (inherits(cholm, "try-error")) return(lik)

  z <- datos - par[1]
  detvarcov <- sum(log(diag(cholm)))

  lik <- 0.5 * (2 * detvarcov + sum(z * backsolve(cholm,
            forwardsolve(cholm, z, transpose = TRUE, upper.tri = TRUE))))

  return(lik)
}

## returns the estimates, the standard errors and the value
fit_mle <- function(datos, D, nu, start = c(0, 0, -3, log(median(D[D > 0])))) {
  o   <- optim(start, loglik, datos = datos, D = D, nu = nu,
               hessian = TRUE, control = list(maxit = 5000, reltol = 1e-10))
  est <- c(mu = o$par[1], sigma2 = exp(o$par[2]),
           nugget = exp(o$par[3]), range = exp(o$par[4]))
  sel <- sqrt(diag(solve(o$hessian)))          # on the working scale
  se  <- sel * c(1, est[2], est[3], est[4])    # delta method
  names(se) <- names(est)
  list(est = est, se = se, loglik = -o$value, conv = o$convergence)
}
