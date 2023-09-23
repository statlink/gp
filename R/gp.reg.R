gp.reg <- function(y, x, tol = 1e-7) {

  x <- model.matrix( y ~ ., data.frame(x) )

  dm <- dim(x)
  n <- dm[1]  ;  k <- dm[2]
  Mod <- Rfast::qpois.reg(x[, -1], y, tol = tol)
  be <- Mod$be
  est <- as.vector( x %*% be )
  mi <- exp( est )
  con <- sum( lgamma(y + 1) )
  pois.loglik <-  - sum(mi) + sum(y * est) - con
  lik1 <- pois.loglik
  f1 <- sqrt( Mod$phi )
  sy <- sum(y)
  y1 <- y - 1
  sx <- Rfast::colsums(x)

  com <- mi + (f1 - 1) * y
  der <- sx + Rfast::eachcol.apply(x, y1 * mi / com ) - Rfast::eachcol.apply(x, mi) / f1
  der2 <- crossprod(x, x * y1 * mi * (f1 - 1) * y / com^2) - 1/f1 * crossprod(x, mi * x)
  be <- be - solve(der2, der)
  est <- as.vector( x %*% be )
  mi <- exp( est )
  lik2 <- sum(est) + sum( y1 * log(com) ) - log(f1) * sy - sum(com) / f1

  while ( lik2 - lik1 > tol ) {
    lik1 <- lik2
    com <- mi + (f1 - 1) * y
    der <- sx + Rfast::eachcol.apply(x, y1 * mi / com ) - Rfast::eachcol.apply(x, mi) / f1
    der2 <- crossprod(x, x * y1 * (f1 - 1) * y * mi / com^2) - 1/f1 * crossprod(x, mi * x)
    be <- be - solve(der2, der)
    est <- as.vector( x %*% be )
    mi <- exp( est )
    lik2 <- sum(est) + sum( y1 * log(com) ) - log(f1) * sy - sum(com) / f1
  }
  f2 <- sqrt( sum( (y - mi)^2/mi ) / (n - k) )

  while ( abs(f1 - f2) > tol ) {
    f1 <- f2
    com <- mi + (f1 - 1) * y
    der <- sx + Rfast::eachcol.apply(x, y1 * mi / com ) - Rfast::eachcol.apply(x, mi) / f1
    der2 <- crossprod(x, x * y1 * mi * (f1 - 1) * y / com^2) - 1/f1 * crossprod(x, mi * x)
    be <- be - solve(der2, der)
    est <- as.vector( x %*% be )
    mi <- exp( est )
    lik2 <- sum(est) + sum( y1 * log(com) ) - log(f1) * sy - sum(com) / f1

    while ( lik2 - lik1 > tol ) {
      lik1 <- lik2
      com <- mi + (f1 - 1) * y
      der <- sx + Rfast::eachcol.apply(x, y1 * mi / com ) - Rfast::eachcol.apply(x, mi) / f1
      der2 <- crossprod(x, x * y1 * (f1 - 1) * y * mi / com^2) - 1/f1 * crossprod(x, mi * x)
      be <- be - solve(der2, der)
      est <- as.vector( x %*% be )
      mi <- exp( est )
      lik2 <- sum(est) + sum( y1 * log(com) ) - log(f1) * sy - sum(com) / f1
    }
    f2 <- sqrt( sum( (y - mi)^2/mi ) / (n - k) )
  }

  list(pois.loglik = pois.loglik, gp.loglik = lik2 - con, be = be, phi = f2)
}











