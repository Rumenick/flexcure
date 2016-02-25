# Marshall-Olkin extended exponential distribution (Marshall & Olkin, 1997)
dmoee <- function(x, mu = 0, alpha = 1, log=FALSE){
  dens <- suppressWarnings(dmoeev(x = x, mu = mu, sigma = 1, alpha = alpha, log = log))
  dens
}

pmoee <- function(q, mu = 0, alpha = 1, lower.tail = TRUE, log.p = FALSE){
  prob <- suppressWarnings(pmoeev(q = q, mu = mu, sigma = 1, alpha = alpha, lower.tail = lower.tail, log.p = log.p))
}

# Extreme value distribuition
dev <- function(x, mu = 0, sigma = 1, log=FALSE){
  dens <- suppressWarnings(dmoeev(x = x, mu = mu, sigma = sigma, alpha = 1, log = log))
  dens
}

pev <- function(q, mu = 0, sigma, lower.tail = TRUE, log.p = FALSE){
  prob <- suppressWarnings(pmoeev(q = q, mu = mu, sigma = sigma, alpha = 1, lower.tail = lower.tail, log.p = log.p))
}

# Likelihood ratio test

LRT <- function(fitg, fits){
  distname <- switch(fits$timedist,
                     exp = "Exponential",
                     weibull = "Weibull",
                     ev = "Extreme value (or Weibull)",
                     lnorm = "Lognormal",
                     gamma = "Gamma",
                     llogis = "Log-logistic",
                     moee = "Marshall-Olkin extended extreme value (or exponential)",
                     moeev = "Marshall-Olkin extended extreme value (or Weibull)",
                     gengamma = "Extended generalized gamma",
                     genf = "Generalized F",
                     "Unknown")
  gl <- attr(logLik(fitg), "df") - attr(logLik(fits), "df")
  est <- 2 * (logLik(fitg) - logLik(fits))
  p <- pchisq(est, gl, lower.tail = F)
  cat("      Likelihood ratio test\n\n")
  pvalue <- if(p < .Machine$double.eps) paste("p-value <", format(.Machine$double.eps, digits = 2)) else paste("p-value =", p)
  cat(" LRS =", est, ",", "df =", gl, ",", pvalue, "\n alternative hypothesis:\n",
      paste(distname, ifelse(is.null(fits$dcure), "", fits$dcure), sep=" "), "\n is not suitable\n")
}

# Plot method for survfit objects using ggplot2

createSurvivalFrame <- function(f.survfit){
  # initialise frame variable
  f.frame <- NULL
  # check if more then one strata
  if(length(names(f.survfit$strata)) == 0){
    # create data.frame with data from survfit
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower)
    # create first two rows (start at 1)
    f.start <- data.frame(time=c(0, f.frame$time[1]), n.risk=c(f.survfit$n, f.survfit$n), n.event=c(0,0),
                          n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1))
    # add first row to dataset
    f.frame <- rbind(f.start, f.frame)
    # remove temporary data
    rm(f.start)
  }
  else {
    # create vector for strata identification
    f.strata <- NULL
    for(f.i in 1:length(f.survfit$strata)){
      # add vector for one strata according to number of rows of strata
      f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
    }
    # create data.frame with data from survfit (create column for strata)
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower, strata=factor(f.strata))
    # remove temporary data
    rm(f.strata)
    # create first two rows (start at 1) for each strata
    for(f.i in 1:length(f.survfit$strata)){
      # take only subset for this strata from data
      f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])
      # create first two rows (time: 0, time of first event)
      f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survfit[f.i]$n, 2), n.event=c(0,0),
                            n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), strata=rep(names(f.survfit$strata)[f.i],
                                                                                                 2))
      # add first two rows to dataset
      f.frame <- rbind(f.start, f.frame)
      # remove temporary data
      rm(f.start, f.subset)
    }
    # reorder data
    f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
    # rename row.names
    rownames(f.frame) <- NULL
  }
  # return frame
  return(f.frame)
}

qplot_survival <- function(f.survfit, f.CI="default", f.shape=3){
  f.frame <- createSurvivalFrame(f.survfit)
  # use different plotting commands dependig whether or not strata's are given
  if("strata" %in% names(f.frame) == FALSE){
    # confidence intervals are drawn if not specified otherwise
    if(f.CI=="default" | f.CI==TRUE){
      # create plot with 4 layers (first 3 layers only events, last layer only censored)
      # hint: censoring data for multiple censoring events at timepoint are overplotted
      # (unlike in plot.survfit in survival package)
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time,
                                                                                            y=upper), direction="hv", linetype=3) + geom_step(aes(x=time,y=lower), direction="hv", linetype=3) +
        geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv), shape=f.shape)
    }
    else {
      # create plot without confidence intervalls
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") +
        geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv), shape=f.shape)
    }
  }
  else {
    if(f.CI=="default" | f.CI==FALSE){
      # without CI
      ggplot(data=f.frame, aes(group=strata, colour=strata)) + geom_step(aes(x=time, y=surv),
                                                                         direction="hv") + geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv), shape=f.shape)
    }
    else {
      # with CI (hint: use alpha for CI)
      ggplot(data=f.frame, aes(colour=strata, group=strata)) + geom_step(aes(x=time, y=surv),
                                                                         direction="hv") + geom_step(aes(x=time, y=upper), direction="hv", linetype=3, alpha=0.5) +
        geom_step(aes(x=time,y=lower), direction="hv", linetype=3, alpha=0.5) +
        geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv), shape=f.shape)
    }
  }
}

