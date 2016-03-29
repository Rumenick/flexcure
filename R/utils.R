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
#' @title Likelihood ratio test
#'
#' @name LRT
#'
#' @description Computes the likelihood ratio test.
#'
#' @param fitg an object that stores the results of
#' curereg fit of the model under the null hypothesis.
#' @param fits an object that stores the results of curereg fit
#' of the model under the alternative hypothesis.
#'
#' @details The objects fitg and fits are obtained
#' using the usual options passed to the curereg
#' function.
#'
#' @return
#' \item{LRS }{
#' the value of the likelihood ratio statistic.
#' }
#' \item{pvalue }{
#' the p value of test under null hypothesis chi-square distribution.
#' }
#'
#' @author Rumenick Pereira da Silva \email{rumenickps@gmail.com}
#'
#' @examples
#'
#' data(e1684)
#' fitg <- curereg(Surv(FAILTIME, FAILCENS) ~ 1, cureformula = ~ 1,
#'                  data = e1684, timedist = "moeev", ncausedist = "bernoulli")
#' fits <- curereg(Surv(FAILTIME, FAILCENS) ~ 1, cureformula = ~ 1,
#'                  data = e1684, timedist = "moee", ncausedist = "bernoulli")
#' LRT(fitg, fits)
#'
#' @keywords test


#' @export
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
  invisible(list(LRS = est, pvalue = pvalue))
}

