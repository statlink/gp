\name{gp.reg}
\alias{gp.reg}
\alias{gp.reg2}
\title{
Generalized Poisson regression
}

\description{
Generalized Poisson regression.
}

\usage{
gp.reg(y, x, tol = 1e-7)
gp.reg2(y, x, tol = 1e-7)
}

\arguments{
\item{y}{
The response variable, a vector with non negative integer values.
}
\item{x}{
A data.frame or a matrix with the independent variables.
}
\item{tol}{
The tolerance value to terminate the optimization.
}
}

\details{
The loglikelihood of the generalised Poisson distribution when covariates are present is the following (Consul & Famoye, 1992):
\deqn{
\ell(\beta, \phi)=\sum_{i=1}^n\log(\mu_i) + \sum_{i=1}^n(y_i-1)\log{[\mu_i+(\phi-1)y_i]}-
\log{\phi}\sum_{i=1}^ny_i-\frac{1}{\phi}\sum_{i=1}^n[\mu_i+(\phi-1)y_i]-\sum_{i=1}^n\log{(y_i)},
}
where \eqn{\mu_i=e^{\sum_{j=0}^kX_{ij}\beta_j}}, \eqn{n} denotes the sample size, \eqn{k} is the number of \eqn{\beta} coefficients, and \eqn{\phi > 0}.

Breslow (1984) suggested the (moment) estimation of a dispersion parameter by equating the chi-square statistic to its degrees of freedom. For the generalised Poisson regression model, this leads to
\eqn{\sum_{i=1}^n\frac{(y_i-\mu_i)^2}{\mu_i\phi^2}=n-k} and we solve this for \eqn{\phi}.

According to Consul and Famoye (1992) we begin by fitting a Poisson regression model and obtain initial values for \eqn{\beta_s} and \eqn{\phi}. If \eqn{\hat{\phi} \approx 1}, it implies that the Poisson regression
model is appropriate and no further estimation needs to be done. However, if \eqn{\hat{\phi} \neq 1}, this is used to obtain new values of the estimated \eqn{\beta_s} by maximizing the log-likelihood. This process is iterated until we obtain a stable solution.

The function as seen below returns the log-likelihood of the initial Poisson regression as well. This is useful if one wants to test, via the log-likelihood ratio test as 1 degree of freedom, if the generalized Poisson regression is to be preferred over the Poisson regression.

gp.reg() estimates the \eqn{beta} coefficients using Newton-Raphson, whereas gp.reg2() uses the \code{\link{optim}} function. For some reason these two do not always agree. One might yield higher log-likelihood than the other and this is why I offer both ways.
}

\value{
A list including:
\item{pois.loglik}{
The initial Poisson regression log-likelihood.
}
\item{gp.loglik}{
The generalized Poisson regression log-likelihood.
}
\item{be}{The estimated \eqn{beta} coefficients.
}
\item{phi}{
The estimated \eqn{\phi} parameter.
}

}

\references{
Consul P.C. & Famoye F. (1992). Generalized poisson regression model. Communications in Statistics - Theory and Methods, 21(1): 89--109.

Breslow N. E. (1984). Extra-Poisson variation in log-linear models. Journal of the Royal Statistical Society: Series C (Applied Statistics), 33(1): 38--44.
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
n <- 500
x <- matrix (rnorm(n * 2), nrow = n, ncol = 2)
be <- c(1, 1)
mi <- x[, 1] * be[1] + x[, 2] * be[2] + 1
mi <- exp(mi)
y <- numeric(n)
for (i in 1:n)  y[i] <- rgp(2, mi[i], 0.5, method = "Inversion")[1]
gp.reg(y, x)
gp.reg2(y, x)
}

