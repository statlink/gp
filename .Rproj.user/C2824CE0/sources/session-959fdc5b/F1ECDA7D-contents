pgp <- function(y, theta, lambda) {

  n <- length(y)
  prob <- numeric(n)
  for (i in 1:n) {
    a <- log( theta ) + (0:y[i] - 1) * log( theta + lambda * (0:y[i]) ) - theta -
        ( 0:y[i] ) * lambda - lgamma( 1:(y[i] + 1) )
    prob[i] <- sum( exp(a) )
  }
  prob

}
