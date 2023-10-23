\name{rgp}
\alias{rgp}

\title{
Random values simulation from the generalized Poisson distribution
}

\description{
Random values simulation from the generalized Poisson distribution.
}

\usage{
rgp(n, theta, lambda, method)
}

\arguments{
\item{n}{
The number of random values to generate.
}
\item{theta}{
The value of the \eqn{\theta} parameter.
}
\item{lambda}{
The value of the \eqn{\lambda} parameter.
}
\item{method}{
The simulation method to use. The available options are:
"Inversion", "Branching", "Normal-Approximation", "Build-Up" and "Chop-Down".
}
}

\details{
The values \eqn{\theta} and \eqn{\lambda} affect the method of simulation to use.
See Li et al. (2020) for more information.
}

\value{
A vector with values from the generalized Poisson distribution.
}

\references{
Li H., Demirtas H. & Chen R. (2020). RNGforGPD: An R Package for Generation of Univariate and Multivariate Generalized Poisson Data. R JOURNAL, 12(2): 173--188.

Demirtas H. (2017). On accurate and precise generation of generalized Poisson variates. Communications in Statistics - Simulation and Computation, 46(1): 489--499.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{gp.mle}
}
}

\examples{
y <-  rgp(1000, 10, 0.5, method = "Inversion")
gp.mle(y)
}


