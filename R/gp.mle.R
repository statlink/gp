gp.mle <- function(y) {

  fun <- function(par, y, y1, n, sy) {
    theta <- exp( par[1] ) ;  lambda <- 1 / ( 1 + exp(- par[2]) )
    f <- n * log( theta ) + sum( y1 * log( theta + lambda * y ) ) -
         n * theta - sy * lambda
    -f
  }

  sy <- sum(y)
  n <- length(y)
  my <- sy / n
  s2 <- Rfast::Var(y)
  y1 <- y - 1

  theta <- log( sqrt(my^3 / s2 ) )
  f <- 1 - sqrt(my / s2)
  lambda <- log(f / (1 - f))
  ini <- c( theta, lambda )

  mod <- optim(ini, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000) )
  lik1 <- mod$value
  mod <- optim(mod$par, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000) )
  lik2 <- mod$value
  while ( lik1 - lik2 > 1e-6 ) {
    lik1 <- lik2
    mod <- optim(mod$par, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000) )
    lik2 <- mod$value
  }
  lik <-  - mod$value - sum( lgamma(y + 1) )
  par <- mod$par
  par <- c( exp( par[1] ), 1 / ( 1 + exp(- par[2]) ), lik )
  names(par) <- c("theta", "lambda", "loglik")
  par
}





