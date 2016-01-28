# Marshall-Olkin extended Weibull distribution (Marshall & Olkin, 1997)

# dmoew <- function(x, shape, scale, alpha, log = FALSE) {
#   if (!check.moew(x=x, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(x)))
#   pdf <- numeric(length(x))
#   pdf[x<=0] <- 0
#   xx <- x[x>0]
#   num <- log(alpha) + log(shape) + shape * (log(xx) - log(scale)) - log(xx) - (xx/scale)^shape
#   den <- -2 * log(1 - (1-alpha) * exp(-(xx/scale)^shape))
#   pdf[x>0] <- num + den
#   if(!log) pdf[x>0] <- exp(pdf[x>0])
#   return(pdf)
# }
#
# pmoew <- function(q, shape, scale, alpha, lower.tail = TRUE, log.p = FALSE) {
#   if (!check.moew(x=q, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(q)))
#   q[q<0] <- 0
#   S <- exp(-(q/scale)^shape)
#   cdf <- (alpha * S)/(1-(1-alpha) * S)
#   if(lower.tail) cdf <- 1 - cdf
#   if(log.p) cdf <- log(cdf)
#   return(cdf)
# }
#
# qmoew <- function(p, shape, scale, alpha, lower.tail = TRUE, log.p = FALSE) {
#   if (!check.moew(x=p, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(p)))
#   if (log.p)
#     p <- exp(p)
#   if (!lower.tail)
#     p <- 1 - p
#   if ((min(p) < 0) || (max(p) > 1))
#     stop("Invalid arguments p")
#   a <- shape
#   b <- scale
#   c <- alpha
#   q <-  b * log((1-(1-c)*p)/(1-p))^(1/a)
#   return(q)
# }
#
# rmoew <- function(n, shape, scale, alpha) {
#   if (!check.moew(x=n, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(n)))
#   if (n <= 0) stop("Invalid arguments n")
#   a <- shape
#   b <- scale
#   c <- alpha
#   u <- runif(n)
#   x <-  b * log((1-(1-c)*u)/(1-u))^(1/a)
#   return(x)
# }
#
# check.moew <- function(x, shape, scale, alpha){
#   ret <- TRUE
#   if ((!is.numeric(shape)) || (!is.numeric(scale)) ||
#       (!is.numeric(alpha)) || (!is.numeric(x)))
#     stop("non-numeric argument to mathematical function")
#   if (any(shape <= 0)) {warning("Non-positive parameter \"shape\""); ret <- FALSE}
#   if (any(scale <= 0)) {warning("Non-positive parameter \"scale\""); ret <- FALSE}
#   if (any(alpha <= 0)) {warning("Non-positive parameter \"alpha\""); ret <- FALSE}
#   ret
# }

# Marshall-Olkin extended extreme value distribution (MOEEV)

dmoeev <- function(x, mu = 0, sigma, alpha = 1, log = FALSE) {
  #if (!check.moew(x=x, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(x)))
  shape <- 1/sigma
  scale <- exp(mu)
  pdf <- numeric(length(x))
  pdf[x<=0] <- 0
  xx <- x[x>0]
  num <- log(alpha) + log(shape) + shape * (log(xx) - log(scale)) - log(xx) - (xx/scale)^shape
  den <- -2 * log(1 - (1-alpha) * exp(-(xx/scale)^shape))
  pdf[x>0] <- num + den
  if(!log) pdf[x>0] <- exp(pdf[x>0])
  return(pdf)
}

pmoeev <- function(q, mu = 0, sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE) {
  #if (!check.moew(x=q, shape=shape, scale=scale, alpha=alpha)) return(rep(NaN, length(q)))
  shape <- 1/sigma
  scale <- exp(mu)
  q[q<0] <- 0
  S <- exp(-(q/scale)^shape)
  cdf <- (alpha * S)/(1-(1-alpha) * S)
  if(lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- log(cdf)
  return(cdf)
}

qmoeev <- function(p, mu = 0, sigma, alpha = 1,  lower.tail = TRUE, log.p = FALSE) {
  #if (!check.moeev(x = p, mu = mu, sigma = sigma, alpha = alpha)) { return(rep(NaN, length(p))) }
  #else{
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    if ((min(p) < 0) || (max(p) > 1)) stop("Invalid arguments p")
    a <- 1/sigma
    b <- exp(mu)
    c <- alpha
    q <-  b * log((1-(1-c)*p)/(1-p))^(1/a)
    return(q)
    #}
}

rmoeev <- function(n, mu = 0, sigma, alpha = 1) {
  if (n <= 0) stop("\nInvalid arguments n")
  #if (!check.moeev(x = n, mu = mu, sigma = sigma, alpha = alpha)) { return(rep(NaN, length(x))) }
  #else{
    a <- 1/sigma
    b <- exp(mu)
    c <- alpha
    u <- runif(n)
    x <-  b * log((1-(1-c)*u)/(1-u))^(1/a)
    return(x)#}
}

