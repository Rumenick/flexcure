% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moeev.R
\name{moeev}
\alias{dmoeev}
\alias{moeev}
\alias{pmoeev}
\alias{qmoeev}
\alias{rmoeev}
\title{Marshall-Olkin extended extreme value (or Weibull) distribution}
\usage{
dmoeev(x, mu = 0, sigma = 1, alpha, log = FALSE)
pmoeev(q, mu = 0, sigma = 1, alpha, lower.tail = TRUE, log.p = FALSE)
qmoeev(p, mu = 0, sigma = 1, alpha, lower.tail = TRUE, log.p = FALSE)
rmoeev(n, mu = 0, sigma = 1, alpha)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mu}{vector of location parameters.}

\item{sigma}{vector of scale parameters.}

\item{alpha}{vector of tilt parameters.}

\item{log, log.p}{logical; if TRUE,  probabilities p are given as log(p).}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If length(n) > 1, the length is taken to be the number
required.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P(X \le x)},
otherwise, \eqn{P(X > x)}.}
}
\value{
\code{dmoeev} gives the density, \code{pmoeev} gives the distribution function, \code{qmoeev} gives the quantile
function, and \code{rmoeev} generates random deviates.
}
\description{
Density, distribution function, quantile function and random generation
for the Marshall-Olkin extended extreme value distribution (MOEEV) distribution with
location parameter mu, scale parameter sigma and tilt paramater alpha.
}
\details{
The Marshall-Olkin Extended Extreme Value (MOEEV) distribution has density

\deqn{f\left(x;\mu,\sigma,\alpha\right)=\frac{\alpha\exp\left[w-e^{w}\right]}{\sigma x\left\{ 1
-\bar{\alpha}\exp\left[-e^{w}\right]\right\} ^{2}, x>0}\]}{f(x;\mu,\sigma,\alpha) = \{\alpha
exp[w-exp(w)]\} / \{\sigma x\{1-(1-\alpha)exp[-exp(w)]\}^2\}, x>0}

where \eqn{w=\frac{\log\left(x\right)-\mu}{\sigma}}{w=[log(x)-\mu]/\sigma} and
\eqn{-\infty<\mu<\infty}, \eqn{\sigma>0} and \eqn{\alpha>0} are the location, scale and tilt
parameters, respectively. If \eqn{\alpha = 1}, we obtain the extreme value distribution.

Consider the parameterisations used by \code{\link{dexp}} and \code{\link{dweibull}}.
If \eqn{\mu = -\log(rate)}{\mu = -log(rate)} and \eqn{\sigma = 1}{\sigma = 1}, we obtain the
Marshall-Olkin Extended Exponential (MOEE) distribution. In the case, that
\eqn{\mu = \log(scale)}{\mu = log(scale)} and \eqn{\sigma = \frac{1}{shape}}{\sigma = 1/shape},
we obtain the Marshall-Olkin Extended Weibull (MOEW) distribution. With the above definitions,

if MOEE:
\code{dmoeev(x, mu = -log(rate), sigma = 1, alpha) = dmoee(x, mu, alpha)}

if MOEW:
\code{dmoeev(x, mu = log(scale), sigma = 1/shape, alpha)}

The Marshall-Olkin extended extreme value distribution simplifies to the exponential and
Weibull distributions with the following parameterisations:

\tabular{lcl}{
\code{dmoeev(x, mu, sigma = 1, alpha = 1)} \tab \code{ = } \tab \code{dexp(x, rate = 1/exp(mu))} \cr
\code{dmoeev(x, mu, sigma, alpha = 1)} \tab \code{ = } \tab  \code{dweibull(x, shape=1/sigma, scale=exp(mu))} \cr
 }
}
\examples{

x <- rmoeev(1000, mu = 2, sigma = 1, alpha = 1)
all.equal(dmoeev(x, mu = 2, sigma = 1, alpha = 1), dexp(x, rate = 1/exp(2)))
x <- rmoeev(1000, mu = 2, sigma = 2, alpha = 1)
all.equal(dmoeev(x, mu = 2, sigma = 2, alpha = 1), dweibull(x, shape=1/2, scale = exp(2)))

}
\author{
Rumenick Pereira da Silva \email{rumenickps@gmail.com}
}
\references{
Marshall, A. W., Olkin, I. (1997). A new method for adding a parameter to a family of
 distributions with application to the Weibull and Weibull families. Biometrika,84(3):641-652.

Marshall, A. W., Olkin, I.(2007). Life Distributions: Structure of Nonparametric,
Semiparametric, and Parametric Families. Springer, New York.
}
\keyword{distribution}

