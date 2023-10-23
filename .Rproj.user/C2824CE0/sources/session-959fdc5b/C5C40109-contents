\name{Cumulative probability mass of the generalized Poisson distribution}
\alias{pgp}
\title{
Cumulative probability mass of the generalized Poisson distribution
}

\description{
Cumulative probability mass of the generalized Poisson distribution.
}

\usage{
pgp(y, theta, lambda)
}

\arguments{
\item{y}{
A vector with non negative integer values.
}
\item{theta}{
The value of the \eqn{\theta} parameter.
}
\item{lambda}{
The value of the \eqn{\lambda} parameter.
}
}

\details{
The cumulative probability mass of the generealized Poisson distribution is computed.
}

\value{
A vector with the probabilities.
}

\references{
Nikoloulopoulos A.K. & Karlis D. (2008). On modeling count data: a comparison of some well-known discrete distributions. Journal of Statistical Computation and Simulation, 78(3): 437--457.

Demirtas H. (2017). On accurate and precise generation of generalized Poisson variates. Communications in Statistics - Simulation and Computation, 46(1): 489--499.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{\link{dgp}, \link{rgp} }
}

\examples{
y <-  rgp(1000, 10, 0.5, method = "Inversion")
a <- gp.mle(y)
pgp(y[1:10], a[1], a[2])
}

