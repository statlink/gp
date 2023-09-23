gp.reg2 <- function(y, x, tol = 1e-7) {

  fun <- function(be, phi, sy, y1) {
    est <- x %*% be
    mi <- exp( est )
    com <- mi + (phi - 1) * y
    f <- sum(est) + sum( y1 * log(com) ) - log(phi) * sy - sum(com) / phi
    -f
  }

  x <- model.matrix( y ~ ., data.frame(x) )

  dm <- dim(x)
  n <- dm[1]  ;  k <- dm[2]
  Mod <- Rfast::qpois.reg(x[, -1], y, tol = tol)
  est <- x %*% Mod$be
  mi <- exp( est )
  con <- sum( lgamma(y + 1) )
  pois.loglik <-  - sum(mi) + sum(y * est) - con
  f1 <- sqrt( Mod$phi )
  sy <- sum(y)
  y1 <- y - 1

  mod <- optim(Mod$be, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
  lik1 <- mod$value
  mod <- optim(mod$par, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
  lik2 <- mod$value
  while ( lik1 - lik2 > tol ) {
    lik1 <- lik2
    mod <- optim(mod$par, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
    lik2 <- mod$value
  }
  mi <- exp( x %*% mod$par )
  f2 <- sqrt( sum( (y - mi)^2/mi ) / (n - k) )

  while ( abs (f1 - f2) > tol ) {
    f1 <- f2
    mod <- optim(mod$par, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
    lik1 <- mod$value
    mod <- optim(mod$par, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
    lik2 <- mod$value
    while ( lik1 - lik2 > tol ) {
      lik1 <- lik2
      mod <- optim(mod$par, fun, phi = f1, sy = sy, y1 = y1, control = list(maxit = 10000) )
      lik2 <- mod$value
    }
    mi <- exp( x %*% mod$par )
    f2 <- sqrt( sum( (y - mi)^2/mi ) / (n - k) )
  }
  be <- mod$par
  gp.loglik <-  - mod$value - con
  list(pois.loglik = pois.loglik, gp.loglik = gp.loglik, be = be, phi = f2)

}






