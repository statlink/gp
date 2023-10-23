\name{gp.mle}
\alias{gp.mle}

\title{
Maximum likelihood estimation of the generalized Poisson distribution
}

\description{
Maximum likelihood estimation of the generalized Poisson distribution.
}

\usage{
gp.mle(y)
}

\arguments{
\item{y}{
A vector with non negative integer values.
}
}

\details{
The probability density function of the generalized Poisson distribution is the following (Nikoloulopoulos & Karlis, 2008):
\deqn{
P(Y=y|\theta, \lambda)=\theta(\theta+\lambda y)^{y-1}\frac{e^{-\theta-\lambda y}}{y!},
\ \ y=0,1... \ \ \theta >0, \ \ 0 \leq \lambda \leq 1.
}
To ensure that \eqn{\theta} is positive we use the "log" link and for \eqn{\lambda} to lie within 0 and 1
we use the "logit" link within the \code{\link{optim}} function.
}

\value{
A vector with three numbers, the \eqn{\theta} and \eqn{\lambda} parameters and the value of the log-likelihood.
}

\references{
Nikoloulopoulos A.K. & Karlis D. (2008). On modeling count data: a comparison of some well-known discrete distributions. Journal of Statistical Computation and Simulation, 78(3): 437--457.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{gp.reg}, \link{rgp}
}
}

\examples{
y <-  rgp(1000, 10, 0.5, method = "Inversion")
gp.mle(y)
}


