## Distribuitons with cure fraction
# dname: Density function
# pname: Cumulative distribution function
# qname: Quantile function (not considered)
# rname: Random number generator


#' @include dpr_generic.R
#' @import stats

## Standard mixture models
# Exponential standard mixture model:

#' @export
dexpsm <- function(x, rate = 1, theta = 0, log=FALSE){
  ret <- log1p(-theta) + dexp(x=x, rate=rate, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pexpsm <- function(q, rate = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pexp(q=q, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rexpsm <- function(n = 1e+3, rate = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rexp',
           rncausedist = 'rbinom',
           partimedist = list(rate = rate),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# Weibull standard mixture model:

#' @export
dweibullsm <- function(x, shape, scale = 1, theta = 0, log=FALSE){
  ret <- log1p(-theta) + dweibull(x=x, shape=shape, scale=scale, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pweibullsm <- function(q, shape, scale = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pweibull(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rweibullsm <- function(n = 1e+3, shape, scale = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rweibull',
           rncausedist = 'rbinom',
           partimedist = list(shape = shape, scale = scale),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# Gamma standard mixture model:

# scale = 1/rate
#' @export
dgammasm <- function(x, shape, rate = 1, scale = 1/rate, theta = 0, log = FALSE){
  ret <- log1p(-theta) + dgamma(x=x, shape=shape, rate=rate, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pgammasm <- function(q, shape, rate = 1, scale = 1/rate, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pgamma(q=q, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rgammasm <- function(n = 1e+3, shape, rate = 1, scale = 1/rate, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgamma',
           rncausedist = 'rbinom',
           partimedist = list(shape = shape, rate = rate, scale = scale),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# log-normal standard mixture model:

#' @export
dlnormsm <- function(x, meanlog = 0, sdlog = 1, theta = 0, log = FALSE){
  ret <- log1p(-theta) + dlnorm(x=x, meanlog = meanlog, sdlog = sdlog, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
plnormsm <- function(q, meanlog = 0, sdlog = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * plnorm(q=q, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rlnormsm <- function(n = 1e+3, meanlog = 0, sdlog = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rlnorm',
           rncausedist = 'rbinom',
           partimedist = list(meanlog = meanlog, sdlog = sdlog),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# Log-logistic standard mixture model:
#' @importFrom flexsurv dllogis
#' @export
dllogissm <- function(x, shape = 1, scale = 1, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dllogis(x=x, shape=shape, scale=scale, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @importFrom flexsurv pllogis
#' @export
pllogissm <- function(q, shape = 1, scale = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * flexsurv::pllogis(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @importFrom flexsurv rllogis
#' @export
rllogissm <- function(n = 1e+3, shape = 1, scale = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rllogis',
           rncausedist = 'rbinom',
           partimedist = list(shape = shape, scale = scale),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# Generalized gamma standard mixture model:

#' @importFrom flexsurv dgengamma
#' @export
dgengammasm <- function(x, mu = 0, sigma = 1, Q, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dgengamma(x = x, mu = mu, sigma = sigma, Q = Q, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @importFrom flexsurv pgengamma
#' @export
pgengammasm <- function(q, mu = 0, sigma = 1, Q, theta = 0, lower.tail = TRUE, log.p=FALSE){
  ret <- (1-theta) * flexsurv::pgengamma(q = q, mu = mu, sigma = sigma, Q = Q, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @importFrom flexsurv rgengamma
#' @export
rgengammasm <- function(n = 1e+3, mu = 0, sigma = 1, Q, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgengamma',
           rncausedist = 'rbinom',
           partimedist = list(mu = mu, sigma = sigma, Q = Q),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# Generalized F standard mixture model:

#' @importFrom flexsurv dgenf
#' @export
dgenfsm <- function(x, mu = 0, sigma = 1, Q, P, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @importFrom flexsurv pgenf
#' @export
pgenfsm <- function(q, mu = 0, sigma = 1, Q, P, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * flexsurv::pgenf(q=q, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @importFrom flexsurv rgenf
#' @export
rgenfsm <- function(n = 1e+3, mu = 0, sigma = 1, Q, P, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgenf',
           rncausedist = 'rbinom',
           partimedist = list(mu = mu, sigma = sigma, Q = Q, P = P),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# MOEEV standard mixture model:

#' @export
dmoeevsm <- function(x, mu, sigma, alpha, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = sigma, alpha = alpha, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pmoeevsm <- function(q, mu, sigma, alpha, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rmoeevsm <- function(n = 1e+3, mu, sigma, alpha, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rmoeev',
           rncausedist = 'rbinom',
           partimedist = list(mu = mu, sigma = sigma, alpha = alpha),
           parncausedist = list(size = 1, prob = 1-theta),
           control = control)
}

# MOEE standard mixture model:

#' @export
dmoeesm <- function(x, mu, alpha, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = 1, alpha = alpha, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pmoeesm <- function(q, mu, alpha, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

# EV standard mixture model:

#' @export
devsm <- function(x, mu, sigma, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = sigma, alpha = 1, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pevsm <- function(q, mu, sigma, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

# Promotion time models
# Exponential promotion time model:

#' @export
dexppt <- function(x, rate = 1, theta = 1, log=FALSE){
  ret <- log(theta) + dexp(x=x, rate=rate, log = TRUE) -
    theta * pexp(q=x, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

#' @export
pexppt <- function(q, rate = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pexp(q=q, rate=rate, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rexppt <- function(n = 1e+3, rate = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rexp',
           rncausedist = 'rpois',
           partimedist = list(rate = rate),
           parncausedist = list(lambda = theta),
           control = control)
}

# Weibull promotion time model:

#' @export
dweibullpt <- function(x, shape, scale = 1, theta = 1, log=FALSE){
  ret <- log(theta) + dweibull(x=x, shape=shape, scale=scale, log = TRUE) -
    theta * pweibull(q=x, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pweibullpt <- function(q, shape, scale = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pweibull(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rweibullpt <- function(n = 1e+3, shape, scale = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rweibull',
           rncausedist = 'rpois',
           partimedist = list(shape = shape, scale = scale),
           parncausedist = list(lambda = theta),
           control = control)
}

# Gamma promotion time model:

# scale = 1/rate
#' @export
dgammapt <- function(x, shape, rate = 1, theta = 1, log = FALSE){
  ret <- log(theta) + dgamma(x=x, shape=shape, rate=rate, log = TRUE) -
    theta * pgamma(q=x, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pgammapt <- function(q, shape, rate = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pgamma(q=q, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rgammapt <- function(n = 1e+3, shape, rate = 1, scale = 1/rate, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgamma',
           rncausedist = 'rpois',
           partimedist = list(shape = shape, rate = rate, scale = scale),
           parncausedist = list(lambda = theta),
           control = control)
}

# Log-normal promotion time model:

#' @export
dlnormpt <- function(x, meanlog = 0, sdlog = 1, theta = 1, log = FALSE){
  ret <- log(theta) + dlnorm(x=x, meanlog = meanlog, sdlog = sdlog, log = TRUE) -
    theta * plnorm(q=x, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
plnormpt <- function(q, meanlog = 0, sdlog = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * plnorm(q=q, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rlnormpt <- function(n = 1e+3, meanlog = 0, sdlog = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rlnorm',
           rncausedist = 'rpois',
           partimedist = list(meanlog = meanlog, sdlog = sdlog),
           parncausedist = list(lambda = theta),
           control = control)
}

# Log-logistic promotion time model:

#' @export
dllogispt <- function(x, shape = 1, scale = 1, theta = 1, log = FALSE){
  ret <- log(theta) + flexsurv::dllogis(x=x, shape=shape, scale=scale, log = TRUE) -
    theta * flexsurv::pllogis(q=x, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pllogispt <- function(q, shape = 1, scale = 1, theta = 1 , lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * flexsurv::pllogis(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rllogispt <- function(n = 1e+3, shape = 1, scale = 1, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rllogis',
           rncausedist = 'rpois',
           partimedist = list(shape = shape, scale = scale),
           parncausedist = list(lambda = theta),
           control = control)
}


# Generalized gamma promotion time model:

#' @export
dgengammapt <- function(x, mu = 0, sigma = 1, Q, theta = 1 , log = FALSE){
  ret <- log(theta) + flexsurv::dgengamma(x=x, mu=mu, sigma=sigma, Q=Q, log = TRUE) -
    theta * flexsurv::pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail = TRUE, log = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pgengammapt <- function(q, mu = 0, sigma = 1, Q, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * flexsurv::pgengamma(q=q, mu=mu, sigma=sigma, Q=Q, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rgengammapt <- function(n = 1e+3, mu = 0, sigma = 1, Q, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgengamma',
           rncausedist = 'rpois',
           partimedist = list(mu = mu, sigma = sigma, Q = Q),
           parncausedist = list(lambda = theta),
           control = control)
}

# Generalized F promotion time model:

#' @export
dgenfpt <- function(x, mu = 0, sigma = 1, Q, P, theta = 1, log = FALSE){
  ret <- log(theta) + flexsurv::dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P, log = TRUE) -
    theta * flexsurv::pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

#' @export
pgenfpt <- function(q, mu = 0, sigma = 1, Q, P, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta*flexsurv::pgenf(q=q, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rgenfpt <- function(n = 1e+3, mu = 0, sigma = 1, Q, P, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rgenf',
           rncausedist = 'rpois',
           partimedist = list(mu = mu, sigma = sigma, Q = Q, P = P),
           parncausedist = list(lambda = theta),
           control = control)
}

# MOEEV promotion time model:

#' @export
dmoeevpt <- function(x, theta, mu, sigma, alpha, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = sigma, alpha = alpha, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

#' @export
pmoeevpt <- function(q, theta, mu, sigma, alpha, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

#' @export
rmoeevpt <- function(n = 1e+3, mu, sigma, alpha, theta, control=1e+3) {
  rgeneric(n = n,
           rtimedist = 'rmoeev',
           rncausedist = 'rpois',
           partimedist = list(mu = mu, sigma = sigma, alpha = alpha),
           parncausedist = list(lambda = theta),
           control = control)
}

# MOEE promotion time model:

#' @export
dmoeept <- function(x, theta, mu, alpha, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = 1, alpha = alpha, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

#' @export
pmoeept <- function(q, theta, mu, alpha, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

# EV promotion time model:

#' @export
devpt <- function(x, theta, mu, sigma, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = sigma, alpha = 1, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

#' @export
pevpt <- function(q, theta, mu, sigma, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}
