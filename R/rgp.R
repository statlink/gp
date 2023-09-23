rgp <- function(n, theta, lambda, method) {
  RNGforGPD::GenUniGpois(theta, lambda, n, details = FALSE, method)$data
}
