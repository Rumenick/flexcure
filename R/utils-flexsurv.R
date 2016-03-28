# Standardised procedure for defining density, cumulative
# distribution and quantile functions for time-to-event
# distributions.

# Christopher Jackson (2015). flexsurv: Flexible Parametric
# Survival and Multi-State Models. R package version 0.7.
# https://CRAN.R-project.org/package=flexsurv.

dbase <- function(dname, lower.tail=TRUE, log=FALSE, ...){
  args <- list(...)
  ## Vectorise all arguments, replicating to length of longest argument
  n <- max(sapply(args, length))
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.fn <- paste("check.",dname,sep="")
  check.ret <- do.call(check.fn, args[-1])
  ret[!check.ret] <- NaN
  for (i in seq_along(args))
    ret[is.nan(args[[i]])] <- NaN
  ## name of first arg is x for PDF, haz, or cum haz, q for CDF and p for quantile function
  stopifnot( !(names(args)[1]=="x" && lower.tail==FALSE))
  if (names(args)[1] %in% c("x","q")){
    x <- args[[1]]
    ## PDF, CDF, hazard and cumulative hazard is 0 for any negative time
    ret[!is.nan(ret) & (x<0)] <- if (lower.tail) { if (log) -Inf else 0 } else { if (log) 0 else 1 }
  }
  if (names(args)[1] == "p") {
    p <- args[[1]]
    if (log) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    args[[1]] <- p
    ret[p < 0 | p > 1] <- NaN
    ## should be 0,Inf for p=0,1, but hopefully always handled anyway
    ## Result is NA if x or a parameter is NA
  }
  ## Result is NA if x or a parameter is NA
  nas <- rep(FALSE, n)
  for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
  ret[nas] <- NA
  ind <- !is.nan(ret) & !nas
  if (names(args)[1] %in% c("x", "q")) ind <- ind & (x>=0)
  ## Any remaining elements of vector are filled in by standard
  ## formula for hazard
  li <- list(ret=ret, ind=ind)
  for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
  c(li, args)
}

### Standardised procedure for defining random sampling functions

rbase <- function(dname, n, ...){
  ## Vectorise all arguments, replicating to sample length
  if (length(n) > 1) n <- length(n)
  args <- list(...)
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.fn <- paste("check.",dname,sep="")
  check.ret <- do.call(check.fn, args)
  ret[!check.ret] <- NaN
  for (i in seq_along(args))
    ret[is.nan(args[[i]])] <- NaN
  nas <- rep(FALSE, n)
  for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
  ret[nas] <- NA
  ind <- !is.nan(ret) & !nas
  li <- list(ret=ret, ind=ind)
  for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
  c(li, args)
}
