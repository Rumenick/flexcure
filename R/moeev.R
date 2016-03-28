# Marshall-Olkin Extended Extreme Value (or Weibull) distribution - MOEEV (Marshall & Olkin, 1997)

#' @title Marshall-Olkin extended extreme value (or Weibull) distribution
#'
#' @name moeev
#'
#' @aliases dmoeev pmoeev qmoeev rmoeev
#'
#' @usage
#' dmoeev(x, mu = 0, sigma, alpha = 1, log = FALSE)
#' pmoeev(q, mu = 0, sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE)
#' qmoeev(p, mu = 0, sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE)
#' rmoeev(n, mu = 0, sigma, alpha = 1)
#'
#' @description Density, distribution function, quantile function and random generation
#' for the Marshall-Olkin extended extreme value distribution (MOEEV) distribution with
#' location parameter mu, scale parameter sigma and tilt paramater alpha.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number
#' required.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param alpha vector of tilt parameters.
#' @param log,log.p logical; if TRUE,  probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le x)},
#' otherwise, \eqn{P(X > x)}.
#'
#'
#' @details The Marshall-Olkin Extended Extreme Value (MOEEV) distribution has density
#'
#' \deqn{f\left(x;\mu,\sigma,\alpha\right)=\frac{\alpha\exp\left[w-e^{w}\right]}{\sigma x\left\{ 1
#' -\bar{\alpha}\exp\left[-e^{w}\right]\right\} ^{2}, x>0}\]}{f(x;\mu,\sigma,\alpha) = \{\alpha
#' exp[w-exp(w)]\} / \{\sigma x\{1-(1-\alpha)exp[-exp(w)]\}^2\}, x>0}
#'
#' where \eqn{w=\frac{\log\left(x\right)-\mu}{\sigma}}{w=[log(x)-\mu]/\sigma} and
#' \eqn{-\infty<\mu<\infty}, \eqn{\sigma>0} and \eqn{\alpha>0} are the location, scale and tilt
#' parameters, respectively. If \eqn{\alpha = 1}, we obtain the extreme value distribution.
#'
#' Consider the parameterisations used by \code{\link{dexp}} and \code{\link{dweibull}}.
#' If \eqn{\mu = -\log(rate)}{\mu = -log(rate)} and \eqn{\sigma = 1}{\sigma = 1}, we obtain the
#' Marshall-Olkin Extended Exponential (MOEE) distribution. In the case, that
#' \eqn{\mu = \log(scale)}{\mu = log(scale)} and \eqn{\sigma = \frac{1}{shape}}{\sigma = 1/shape},
#' we obtain the Marshall-Olkin Extended Weibull (MOEW) distribution. With the above definitions,
#'
#' if MOEE:
#' \code{dmoeev(x, mu = -log(rate), sigma = 1, alpha)}
#'
#' if MOEW:
#' \code{dmoeev(x, mu = log(scale), sigma = 1/shape, alpha)}
#'
#' The Marshall-Olkin extended extreme value distribution simplifies to the exponential and
#' Weibull distributions with the following parameterisations:
#'
#' \tabular{lcl}{
#' \code{dmoeev(x, mu, sigma = 1, alpha = 1)} \tab \code{ = } \tab \code{dexp(x, rate = 1/exp(mu))} \cr
#' \code{dmoeev(x, mu, sigma, alpha = 1)} \tab \code{ = } \tab  \code{dweibull(x, shape=1/sigma, scale=exp(mu))} \cr
#'  }
#'
#' @return \code{dmoeev} gives the density, \code{pmoeev} gives the distribution function, \code{qmoeev} gives the quantile
#' function, and \code{rmoeev} generates random deviates.
#'
#' @references Marshall, A. W., Olkin, I. (1997). A new method for adding a parameter to a family of
#'  distributions with application to the Weibull and Weibull families. Biometrika,84(3):641-652.
#'
#' Marshall, A. W., Olkin, I.(2007). Life Distributions: Structure of Nonparametric,
#' Semiparametric, and Parametric Families. Springer, New York.
#'
#' @author Rumenick Pereira da Silva <rumenickbf@hotmail.com>
#'
#' @examples
#'
#' x <- rmoeev(1000, mu = 2, sigma = 1, alpha = 1)
#' all.equal(dmoeev(x, mu = 2, sigma = 1, alpha = 1), dexp(x, rate = 1/exp(2)))
#' x <- rmoeev(1000, mu = 2, sigma = 2, alpha = 1)
#' all.equal(dmoeev(x, mu = 2, sigma = 2, alpha = 1), dweibull(x, shape=1/2, scale = exp(2)))
#'
#' @keywords distribution
#'
#' @include utils-flexsurv.R
#' @export
dmoeev <- function(x, mu = 0, sigma, alpha = 1, log = FALSE) {
  d <- dbase("moeev", log = log, x = x, mu = mu, sigma = sigma, alpha = alpha)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  shape <- 1/sigma
  scale <- exp(mu)
  num <- log(alpha) + log(shape) + shape * (log(x) - log(scale)) - log(x) - (x/scale)^shape
  den <- -2 * log(1 - (1-alpha) * exp(-(x/scale)^shape))
  logdens <- num + den
  ret[ind] <- if (log) logdens else exp(logdens)
  return(ret)
}

#' @export
pmoeev <- function(q, mu = 0, sigma = 1, alpha, lower.tail = TRUE, log.p = FALSE) {
  d <- dbase("moeev", lower.tail = lower.tail, log = log.p, q = q, mu = mu, sigma = sigma, alpha = alpha)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  shape <- 1/sigma
  scale <- exp(mu)
  S <- exp(-(q/scale)^shape)
  prob <- (alpha * S)/(1-(1-alpha) * S)
  if(lower.tail) prob <- 1 - prob
  if(log.p) prob <- log(prob)
  ret[ind] <- prob
  return(ret)
}

#' @export
qmoeev <- function(p, mu = 0, sigma = 1, alpha,  lower.tail = TRUE, log.p = FALSE) {
  d <- dbase("moeev", lower.tail = lower.tail, log = log.p, p = p, mu = mu, sigma = sigma, alpha = alpha)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  shape <- 1/sigma
  scale <- exp(mu)
  ret[ind] <- scale * log((1-(1-alpha)*p)/(1-p))^(1/shape)
  return(ret)
}

#' @export
rmoeev <- function(n, mu = 0, sigma = 1, alpha) {
  d <- rbase("moeev", n = n, mu = mu, sigma = sigma, alpha = alpha)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  shape <- 1/sigma
  scale <- exp(mu)
  u <- runif(length(mu))
  ret[ind] <-  scale * log((1-(1-alpha)*u)/(1-u))^(1/shape)
  return(ret)
}

check.moeev <- function(mu = 0, sigma = 1, alpha){
  ret <- rep(TRUE, length(mu))
  if (missing(alpha)) stop("Tilt parameter \"alpha\" not given")
  if (any(sigma < 0)) { # no warning for sigma 0, since may occur in optimisation.
    warning("Negative scale parameter \"sigma\"")
    ret[sigma < 0] <- FALSE
  }
  if (any(alpha < 0)) { # no warning for alpha 0, since may occur in optimisation.
    warning("Negative till parameter \"alpha\"")
    ret[alpha < 0] <- FALSE
    }
  ret
}
