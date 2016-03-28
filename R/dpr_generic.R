dgeneric <- function() {

}

pgeneric <- function() {

}

rgeneric <- function(n, rtimedist, rncausedist, partimedist, parncausedist, control) {
  if (n <= 0) stop('\nInvalid arguments n')
  partimedist <- lapply(partimedist, rep_len, length.out = n)
  parncausedist <- lapply(parncausedist, rep_len, length.out = n)
  parncausedist$n <- n

  mi <- do.call(rncausedist, parncausedist)

  ti <- ci <- ri <- runif(n, 0, control)
  partimedist <- lapply(partimedist, rep, times = mi)
  partimedist$n <- sum(mi)
  ri[mi>0] <- sapply(split(do.call(rtimedist, partimedist), rep(1:n, mi)), min)

  ti <- pmin(ri, ci)
  di <- as.numeric(ri < ci)
  return(list(time = ti, status = di, ncause = mi))
}
