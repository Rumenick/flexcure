## Distribuitons with cure fraction
# dname: Density function
# pname: Cumulative distribution function
# qname: Quantile function (not considered)
# rname: Random number generator

# Standard mixture models
# Exponential standard mixture model:
compiler::enableJIT(3)
dexpsm <- function(x, rate = 1, theta = 0, log=FALSE){
  ret <- log1p(-theta) + dexp(x=x, rate=rate, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pexpsm <- function(q, rate = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pexp(q=q, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rexpsm <- function(n = 1e+3, rate = 1, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rexp(n=sum(M), rate=rate)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Weibull standard mixture model:

dweibullsm <- function(x, shape, scale = 1, theta = 0, log=FALSE){
  ret <- log1p(-theta) + dweibull(x=x, shape=shape, scale=scale, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pweibullsm <- function(q, shape, scale = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pweibull(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rweibullsm <- function(n = 1e+3, shape, scale = 1, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rweibull(n=sum(M), shape=shape, scale=scale)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Gamma standard mixture model:

# scale = 1/rate
dgammasm <- function(x, shape, rate = 1, theta = 0, log = FALSE){
  ret <- log1p(-theta) + dgamma(x=x, shape=shape, rate=rate, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pgammasm <- function(q, shape, rate = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pgamma(q=q, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgammasm <- function(n = 1e+3, shape, rate = 1, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rgamma(n=sum(M), shape=shape, rate=rate)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# log-normal standard mixture model:

dlnormsm <- function(x, meanlog = 0, sdlog = 1, theta = 0, log = FALSE){
  ret <- log1p(-theta) + dlnorm(x=x, meanlog = meanlog, sdlog = sdlog, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

plnormsm <- function(q, meanlog = 0, sdlog = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * plnorm(q=q, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rlnormsm <- function(n = 1e+3, meanlog = 0, sdlog = 1, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rlnorm(n=sum(M),meanlog = meanlog, sdlog = sdlog)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Log-logistic standard mixture model:

dllogissm <- function(x, shape = 1, scale = 1, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dllogis(x=x, shape=shape, scale=scale, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pllogissm <- function(q, shape = 1, scale = 1, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * flexsurv::pllogis(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rllogissm <- function(n = 1e+3, shape = 1, scale = 1, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- flexsurv::rllogis(n=sum(M), shape=shape, scale=scale)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Generalized gamma standard mixture model:

dgengammasm <- function(x, mu = 0, sigma = 1, Q, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dgengamma(x = x, mu = mu, sigma = sigma, Q = Q, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pgengammasm <- function(q, mu = 0, sigma = 1, Q, theta = 0, lower.tail = TRUE, log.p=FALSE){
  ret <- (1-theta) * flexsurv::pgengamma(q = q, mu = mu, sigma = sigma, Q = Q, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgengammasm <- function(n = 1e+3, mu = 0, sigma = 1, Q, theta, control=1e+3){
  theta <- rep_len(theta, length.out = n)
  mu <- rep_len(mu, length.out = n)
  sigma <- rep_len(sigma, length.out = n)
  Q <- rep(Q, length.out = n)
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- flexsurv::rgengamma(n=sum(M), mu = mu [M==1], sigma = sigma[M==1], Q = Q[M==1])
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Generalized F standard mixture model:

dgenfsm <- function(x, mu = 0, sigma = 1, Q, P, theta = 0, log = FALSE){
  ret <- log1p(-theta) + flexsurv::dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pgenfsm <- function(q, mu = 0, sigma = 1, Q, P, theta = 0, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * flexsurv::pgenf(q=q, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgenfsm <- function(n = 1e+3, mu = 0, sigma = 1, Q, P, theta, control=1e+3){
  theta <- rep_len(theta, length.out = n)
  mu <- rep_len(mu, length.out = n)
  sigma <- rep_len(sigma, length.out = n)
  Q <- rep(Q, length.out = n)
  P <- rep(P, length.out = n)
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- flexsurv::rgenf(n=sum(M), mu = mu[M==1], sigma = sigma[M==1], Q = Q[M==1], P = P[M==1])
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# rgenfsm <- function(n = 1e+3, mu = 0, sigma = 1, s1, s2, theta, control=1e+3){
#   M <- rbinom(n=n, size=1, prob=1-theta)
#   C <- runif(n=n, min=0, max=control)
#   T[M==1] <- rgf(n=sum(M), mu=mu, sigma=sigma, s1=s1, s2=s2)
#   T[M==0] <- C[M==0]
#   y <- apply(cbind(T,C), 1, min)
#   d <- rep(0, n)
#   d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
#   return(list(time = y, status = d, cure = M))
# }

# MOEEV standard mixture model:

dmoeevsm <- function(x, mu, sigma, alpha, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = sigma, alpha = alpha, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pmoeevsm <- function(q, mu, sigma, alpha, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rmoeevsm <- function(n = 1e+3, mu, sigma, alpha, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rmoeev(n=sum(M), mu = mu, sigma = sigma, alpha = alpha)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# MOEE standard mixture model:

dmoeesm <- function(x, mu, alpha, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = 1, alpha = alpha, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pmoeesm <- function(q, mu, alpha, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rmoeesm <- function(n = 1e+3, mu, alpha, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rmoeev(n=sum(M), mu = mu, sigma = 1, alpha = alpha)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# EV standard mixture model:

devsm <- function(x, mu, sigma, theta, log = FALSE){
  ret <- log1p(-theta) + dmoeev(x = x,  mu = mu, sigma = sigma, alpha = 1, log = TRUE)
  if (!log) ret <- exp(ret)
  ret
}

pevsm <- function(q, mu, sigma, theta, lower.tail = TRUE, log.p = FALSE){
  ret <- (1-theta) * pmoeev(q = q,  mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p=FALSE)
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

revsm <- function(n = 1e+3, mu, sigma, theta, control=1e+3){
  M <- rbinom(n=n, size=1, prob=1-theta)
  C <- runif(n=n, min=0, max=control)
  T[M==1] <- rmoeev(n=sum(M), mu = mu, sigma = sigma, alpha = 1)
  T[M==0] <- C[M==0]
  y <- apply(cbind(T,C), 1, min)
  d <- rep(0, n)
  d[M==1] <- ifelse(T[M==1] > C[M==1], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Promotion time models
# Exponential promotion time model:

dexppt <- function(x, rate = 1, theta = 1, log=FALSE){
  ret <- log(theta) + dexp(x=x, rate=rate, log = TRUE) -
    theta * pexp(q=x, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

pexppt <- function(q, rate = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pexp(q=q, rate=rate, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rexppt <- function(n = 1e+3, rate = 1, theta,  control = 1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rexp(n=M, rate=rate)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# Weibull promotion time model:

dweibullpt <- function(x, shape, scale = 1, theta = 1, log=FALSE){
  ret <- log(theta) + dweibull(x=x, shape=shape, scale=scale, log = TRUE) -
    theta * pweibull(q=x, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

pweibullpt <- function(q, shape, scale = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pweibull(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rweibullpt <- function(n = 1e+3, shape = 1, scale = 1, theta = 0.1, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rweibull(n=M, shape=shape, scale=scale)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# Gamma promotion time model:

# scale = 1/rate
dgammapt <- function(x, shape, rate = 1, theta = 1, log = FALSE){
  ret <- log(theta) + dgamma(x=x, shape=shape, rate=rate, log = TRUE) -
    theta * pgamma(q=x, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

pgammapt <- function(q, shape, rate = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * pgamma(q=q, shape=shape, rate=rate, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgammapt <- function(n = 1e+3, shape = 1, rate = 1, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rgamma(n=M, shape=shape, rate=rate)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# Log-normal promotion time model:

dlnormpt <- function(x, meanlog = 0, sdlog = 1, theta = 1, log = FALSE){
  ret <- log(theta) + dlnorm(x=x, meanlog = meanlog, sdlog = sdlog, log = TRUE) -
    theta * plnorm(q=x, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

plnormpt <- function(q, meanlog = 0, sdlog = 1, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * plnorm(q=q, meanlog = meanlog, sdlog = sdlog, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rlnormpt <- function(n = 1e+3, meanlog = 0, sdlog = 1, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rlnorm(n=M,  meanlog = meanlog, sdlog = sdlog)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# Log-logistic promotion time model:

dllogispt <- function(x, shape = 1, scale = 1, theta = 1, log = FALSE){
  ret <- log(theta) + flexsurv::dllogis(x=x, shape=shape, scale=scale, log = TRUE) -
    theta * flexsurv::pllogis(q=x, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

pllogispt <- function(q, shape = 1, scale = 1, theta = 1 , lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * flexsurv::pllogis(q=q, shape=shape, scale=scale, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rllogispt <- function(n = 1e+3, shape = 1, scale = 1, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(flexsurv::rllogis(n=M, shape=shape, scale=scale)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# Generalized gamma promotion time model:

dgengammapt <- function(x, mu = 0, sigma = 1, Q, theta = 1 , log = FALSE){
  ret <- log(theta) + flexsurv::dgengamma(x=x, mu=mu, sigma=sigma, Q=Q, log = TRUE) -
    theta * flexsurv::pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail = TRUE, log = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

pgengammapt <- function(q, mu = 0, sigma = 1, Q, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta * flexsurv::pgengamma(q=q, mu=mu, sigma=sigma, Q=Q, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgengammapt <- function(n = 1e+3, mu = 0, sigma = 1, Q, theta, control=1e+3){
  theta <- rep_len(theta, length.out = n)
  mu <- rep_len(mu, length.out = n)
  sigma <- rep_len(sigma, length.out = n)
  Q <- rep(Q, length.out = n)
  M <- rpois(n = n, lambda = theta)
  I0  <- M!=0
  theta <- rep(theta[I0], M[I0])
  mu <- rep(mu[I0], M[I0])
  sigma <- rep(sigma[I0], M[I0])
  Q <- rep(Q[I0], M[I0])
  t <- numeric(n)
  C <- runif(n=n, min=0, max=control)
  t[I0] <- sapply(X = split(flexsurv::rgengamma(n = sum(M), mu=mu, sigma=sigma, Q=Q),
                            rep((1:n)[I0], M[I0])),
                  FUN = min)
  t[!I0] <- C[!I0]
  y <- apply(cbind(t, C), 1, min)
  d <- rep(0, n)
  d[I0] <- ifelse(t[I0] > C[I0], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# Generalized F promotion time model:

dgenfpt <- function(x, mu = 0, sigma = 1, Q, P, theta = 1, log = FALSE){
  ret <- log(theta) + flexsurv::dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P, log = TRUE) -
    theta * flexsurv::pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE)
  if (!log) ret <- exp(ret)
  ret
}

pgenfpt <- function(q, mu = 0, sigma = 1, Q, P, theta = 1, lower.tail = TRUE, log.p = FALSE){
  ret <- 1 - exp(-theta*flexsurv::pgenf(q=q, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rgenfpt <- function(n = 1e+3, mu = 0, sigma = 1, Q, P, theta = 0.1, control=1e+3){
  theta <- rep_len(theta, length.out = n)
  mu <- rep_len(mu, length.out = n)
  sigma <- rep_len(sigma, length.out = n)
  Q <- rep(Q, length.out = n)
  P <- rep(P, length.out = n)
  M <- rpois(n = n, lambda = theta)
  I0  <- M!=0
  theta <- rep(theta[I0], M[I0])
  mu <- rep(mu[I0], M[I0])
  sigma <- rep(sigma[I0], M[I0])
  Q <- rep(Q[I0], M[I0])
  P <- rep(P[I0], M[I0])
  t <- numeric(n)
  C <- runif(n=n, min=0, max=control)
  t[I0] <- sapply(X = split(flexsurv::rgenf(n = sum(M), mu=mu, sigma=sigma, Q=Q, P=P),
                            rep((1:n)[I0], M[I0])),
                  FUN = min)
  t[!I0] <- C[!I0]
  y <- apply(cbind(t, C), 1, min)
  d <- rep(0, n)
  d[I0] <- ifelse(t[I0] > C[I0], 0, 1)
  return(list(time = y, status = d, cure = M))
}

# MOEEV promotion time model:

dmoeevpt <- function(x, theta, mu, sigma, alpha, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = sigma, alpha = alpha, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

pmoeevpt <- function(q, theta, mu, sigma, alpha, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = sigma, alpha = alpha, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rmoeevpt <- function(n = 1e+3, mu, sigma, alpha, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rmoeev(n=M, mu = mu, sigma = sigma, alpha = alpha)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# MOEE promotion time model:

dmoeept <- function(x, theta, mu, alpha, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = 1, alpha = alpha, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

pmoeept <- function(q, theta, mu, alpha, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = 1, alpha = alpha, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

rmoeept <- function(n = 1e+3, mu, alpha, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rmoeev(n=M, mu = mu, sigma = 1, alpha = alpha)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

# EV promotion time model:

devpt <- function(x, theta, mu, sigma, log = FALSE){
  ret <- log(theta) + dmoeev(x = x, mu = mu, sigma = sigma, alpha = 1, log = TRUE) -
    theta * pmoeev(q = x, mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE)
  if(!log) ret <- exp(ret)
  ret
}

pevpt <- function(q, theta, mu, sigma, lower.tail = TRUE, log.p=FALSE){
  ret <- 1 - exp(-theta * pmoeev(q = q, mu = mu, sigma = sigma, alpha = 1, lower.tail = TRUE, log.p = FALSE))
  if (!lower.tail) ret <- 1-ret
  if (log.p) ret <- log1p(ret-1)
  ret
}

revpt <- function(n = 1e+3, mu, sigma, theta, control=1e+3){
  M <- matrix(rpois(n=n, lambda=theta))
  C <- runif(n=n, min=0, max=control)
  f <- function(M){ ifelse(M > 0, min(rmoeev(n=M, mu = mu, sigma = sigma, alpha = 1)), 0) }
  T <- apply(M, 1, f); T[M==0] <- C[M==0]
  y <- apply(cbind(T, C), 1, min)
  d <- ifelse(T < C, 1, 0)
  return(list(time = y, status = d, cure = M[,1]))
}

