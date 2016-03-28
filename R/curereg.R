## Parameter transformations

# logit function
# if inverse = TRUE
# inverse logit function
logit <- function(x, bvalue=.Machine$double.eps, inverse = FALSE){
  if (!inverse && length(bvalue)){
    x[x <= 0] <- bvalue
    x[x >= 1] <- 1-bvalue
  }
  output <- ifelse(!inverse, expression(log(x) - log1p(-x)), expression(exp(x - log1p(exp(x)))))
  return(eval(output))
}

# log function
# if inverse = TRUE
# exponential function
loge <- function(x, bvalue=.Machine$double.eps, inverse = FALSE){
  if (!inverse && length(bvalue)){ x[x <= 0] <- bvalue}
  output <- ifelse(!inverse,  expression(log(x)), expression(exp(x)))
  return(eval(output))
}

## list of avaliable distribuition:
# Standart Mixture - SM (sm)
# Promotion Time - PT (pt)

curereg.dists <- list(
  expsm = list(
    name = "expsm",
    pars = c("theta", "rate"),
    location = c("rate"),
    transforms = c(logit, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      c(initp, 1/mean(t))
    }}),
  exppt = list(
    name = "exppt",
    pars = c("theta", "rate"),
    location = c("rate"),
    transforms = c(loge, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      c(-log(initp), 1/mean(t))
    }}),
  weibullsm = list(
    name = "weibullsm",
    pars = c("theta", "shape", "scale"),
    location = c("scale"),
    transforms = c(logit, loge, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(initp, 1, exp(mean(lt) + 0.572))
    }}),
  weibullpt = list(
    name = "weibullpt",
    pars = c("theta","shape","scale"),
    location = c("scale"),
    transforms = c(loge, loge, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp), 1, exp(mean(lt) + 0.572))
    }}),
  gammasm = list(
    name = "gammasm",
    pars = c("theta", "shape", "rate"),
    location = c("rate"),
    transforms = c(logit, loge, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      m = mean(t)
      v = var(t)
      c(initp, m^2/v, m/v)
    }}),
  gammapt = list(
    name = "gammapt",
    pars = c("theta","shape","rate"),
    location = c("rate"),
    transforms = c(loge, loge, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      m = mean(t)
      v = var(t)
      c(-log(initp), m^2/v, m/v)
    }}),
  lnormsm = list(
    name = "lnormsm",
    pars = c("theta", "meanlog", "sdlog"),
    location = c("meanlog"),
    transforms = c(logit, identity, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(initp, mean(lt), sd(lt))
    }}),
  lnormpt = list(
    name = "lnormpt",
    pars = c("theta", "meanlog", "sdlog"),
    location = c("meanlog"),
    transforms = c(loge, identity, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp), mean(lt), sd(lt))
    }}),
  llogissm = list(
    name = "llogissm",
    pars = c("theta", "shape", "scale"),
    location = c("scale"),
    transforms = c(logit, loge, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(initp, 1, exp(median(lt)))
    }}),
  llogispt = list(
    name = "llogispt",
    pars = c("theta", "shape", "scale"),
    location = c("scale"),
    transforms = c(loge, loge, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp), 1, exp(median(lt)))
    }}),
  gengammasm = list(
    name = "gengammasm",
    pars = c("theta","mu","sigma","Q"),
    location = c("mu"),
    transforms = c(logit, identity, loge, identity),
    inv.transforms = c(function(x)logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), identity),
    inits = function(initp) { function(t){
      lt <- log(t[t > 0])
      c(initp, mean(lt), sd(lt), 0)
    }}),
  gengammapt = list(
    name = "gengammapt",
    pars = c("theta", "mu", "sigma", "Q"),
    location = c("mu"),
    transforms = c(loge, identity, loge, identity),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), identity),
    inits = function(initp) { function(t){
      lt <- log(t[t > 0])
      c(-log(initp), mean(lt), sd(lt), 0)
    }}),
  genfsm = list(
    name = "genfsm",
    pars = c("theta","mu","sigma","Q","P"),
    location = c("mu"),
    transforms = c(logit, identity, loge, identity, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function(t){
      lt <- log(t[t > 0])
      c(initp, mean(lt), sd(lt), 0, 1)
    }}),
  genfpt = list(
    name = "genfpt",
    pars = c("theta", "mu", "sigma", "Q", "P"),
    location = c("mu"),
    transforms = c(loge, identity, loge, identity, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), identity,function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function(t){
      lt <- log(t[t > 0])
      c(-log(initp), mean(lt), sd(t), 0, 1)
    }}),
  ev = list(
    name = "ev",
    pars = c("mu", "sigma"),
    location = c("mu"),
    transforms = c(identity, loge),
    inv.transforms = c(identity, function(x) loge(x, inverse = TRUE)),
    inits = function (t)
    {
      lt <- log(t[t > 0])
      c(mean(lt) + 0.572, sd(lt))
    }),
  moee = list(
    name = "moee",
    pars = c("mu", "alpha"),
    location = c("mu"),
    transforms = c(identity, loge),
    inv.transforms = c(identity, function(x) loge(x, inverse = TRUE)),
    inits = function (t)
    {
      lt <- log(t[t > 0])
      c(mean(lt) + 0.572, 1)
    }),
  moeev = list(
    name = "moeev",
    pars = c("mu", "sigma", "alpha"),
    location = c("mu"),
    transforms = c(identity, loge, loge),
    inv.transforms = c(identity, function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function (t)
    {
      lt <- log(t[t > 0])
      c(mean(lt) + 0.572, 1, 1)
    }),
  evsm = list(
    name = "evsm",
    pars = c("theta", "mu", "sigma"),
    location = c("mu"),
    transforms = c(logit, identity, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(initp, mean(lt) + 0.572, 1)
    }}),
  moeesm = list(
    name = "moeesm",
    pars = c("theta", "mu", "alpha"),
    location = c("mu"),
    transforms = c(logit, identity, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(initp, mean(lt) + 0.572, 1)
    }}),
  moeevsm = list(
    name = "moeevsm",
    pars = c("theta", "mu", "sigma", "alpha"),
    location = c("mu"),
    transforms = c(logit, identity, loge, loge),
    inv.transforms = c(function(x) logit(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
    lt <- log(t[t > 0])
    c(initp, mean(lt) + 0.572, 1, 1)
    }}),
  evpt = list(
    name = "evpt",
    pars = c("theta", "mu", "sigma"),
    location = c("mu"),
    transforms = c(loge, identity, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp),  mean(lt) + 0.572, 1)
    }}),
  moeept = list(
    name = "moeept",
    pars = c("theta", "mu", "alpha"),
    location = c("mu"),
    transforms = c(loge, identity, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp),  mean(lt) + 0.572, 1)
    }}),
  moeevpt = list(
    name = "moeevpt",
    pars = c("theta", "mu", "sigma", "alpha"),
    location = c("mu"),
    transforms = c(loge, identity, loge, loge),
    inv.transforms = c(function(x) loge(x, inverse = TRUE), identity, function(x) loge(x, inverse = TRUE), function(x) loge(x, inverse = TRUE)),
    inits = function(initp) { function (t)
    {
      lt <- log(t[t > 0])
      c(-log(initp),  mean(lt) + 0.572, 1, 1)
    }})
)

## Function for adjust regression for a flexible parametric survival models with cure fraction:

#' @title Parametric regression models with cure fraction for survival data
#'
#' @name curereg
#'
#' @aliases curereg
#'
#' @description \code{curereg} fits parametric regression models with cure fraction for survival
#' data. This function extends the \code{\link{flexsurvreg}} by the inclusion of the cure fraction
#' in the formulation and adds the Marshall-Olkin extreme value distribution in the comprehensive
#' roll of parametric distributions avaliable.
#'
#' @param formula an object of class "\code{\link{formula}}" which expresses the model to be fitted.
#' The response should be specified as an object of class \code{survival} obtained via
#' \code{\link{Surv}} function and can be set as right, left and interval censored data.
#' The \code{~} separates the response from the covariates which should be specified on the right
#' side, for instance \code{formula = Surv(time, status) ~ age + sex}.
#'
#' @param cureformula a formula defining the cure rate model. In the case of no covariate effect
#' on the cure rate, set \code{cureformula = ~ 1} (default), and for instance
#' \code{cureformula = ~ age + sex} to include coefficients \code{age} and \code{sex} additively.
#' If \code{cureformula = }  or {\code{cureformula = NULL}} the model will be fitted without any
#' cure rate.
#'
#' @param data the data set of class \code{data.frame} or \code{list} which includes all the
#' objects defined in \code{formula} and \code{cureformula}. In case unspecified \code{data},
#' the variables shoud be available in the workspace (see \code{(.GlobalEnv)}).
#'
#' @param weights optional prior weights for the data.
#'
#' @param timedist survival distribution for the non-cured individuals. This can be set as:
#' \code{"exp"} (exponential), \code{"weibull"} (Weibull), \code{"ev"} (extreme value),
#' \code{"gamma"} (Gamma),\code{"lnorm"} (Log-normal), \code{"llogist"} (Log-logistic),
#' \code{"moee"} (extreme value (or exponential) in the Marshall-Olkin family),
#' \code{"moeev"} (extreme value (or Weibull) in the Marshall-Olkin family, the default),
#' \code{"gengamma"} (Generalised Gamma), \code{"genf"} (Generalised F).
#'
#' The exponential, Weibull, log-normal and log-logistic distributions have the same
#' parameterization defined in \code{\link{dexp}}, \code{\link{dweibull}}, \code{\link{dlnorm}}
#' from package \pkg{base} and \code{\link{dllogist}} from package \pkg{flexsurv}, respectively.
#' respectively. These differ from the parametrization used in the package \pkg{survreg}.
#' The generalised Gamma and Generalised F distributions follows the parametrisation in
#' (\code{\link{dgengamma}}) and (\code{\link{dgenf}}), respectively, both available in
#' \pkg{flexsurv}. For the Marshall-Olkin extreme value distribution see \code{\link{dmoeev}}.
#'
#' @param ncausedist distribution of the number of competing causes of the event. This can be
#' set as \code{ncausedist = "bernoulli"} for the standard mixture model and
#' \code{ncausedist = "poisson"} (default), for the promotion time model.
#'
#' @param inits optional list with the initial values for the parameters. This list should be
#' set as \code{inits = list(coef_cure = c(...), coef_time = c(...), sigma = ..., Q = ...)}
#' where \code{coef_cure} is the vector of coefficients for the cure model, \code{coef_time}
#' is the vector of coefficients for the survival regression model, \code{sigma} is  ...
#' and \code{Q} is ... .
#'
#' @param subset optional numeric vector specifying the subset observations from the full data set.
#'
#' @param na.action a function indicating what should happen when \code{NA}'s occur, with
#' possible arguments \code{na.omit} and \code{na.fail}. The default is set by the
#' \code{na.action} setting in \code{options()}.
#'
#' @param method The optimisation method to be used. \code{method = "BFGS"} is the default
#' however "Nelder-Mead", "CG", "L-BFGS-B" and "SANN" can also be used, For more information
#' about the optimisation methods, see \code{\link{optim}}.
#'
#' @param \dots optional arguments for the \code{\link{optim}} and \code{\link{flexsurvreg}}
#' functions.
#'
#' For situations where the default \code{optim} arguments results in lack of converge, consider
#' use \code{control=list(fnscale = value)} with \code{value} a tolerance value with the same
#' magnitude of the log-likelihood function. Usually that happens when the Hessian matrix is
#' not positive definite at some step of the numerical optimisation. An useful tool detect this
#' and verify a possible "slower" convergence is add an appropriate value for \code{trace} to
#' the \code{control}. See \code{\link{optim}} for more information.
#'
#' The argument \code{fixedpars} of the \code{flexsurvreg} function allows the user to input a
#' vector of indices representing the fixed parameters. The arbitrary values for those fixed
#' parameters should have been specified in \code{inits} argument and remain fixed throughout
#' the estimation process.
#'
#' @details Note that the arguments \code{ncausedist} and \code{timedist} set up the model to
#' be fitted. This means that if \code{ncausedist = "poisson"} and \code{timedist = "genf"}
#' the fitted model is obtained considering the promotion-time model with generalised F
#' responses. The improper density function of this model is available in \code{\link{dgenfpt}}.
#' If \code{ncausedist = "bernoulli"} and \code{timedist = "genf"} has been set then the fitted
#' model is calculated considering the standard mixture model with generalised F responses.
#' The improper density function of this model is available in \code{\link{dgenfms}}.
#' \code{d____ms} and \code{p____ms} correspond to the improper density and probability
#' functions for standard mixture models, respectively, where \code{____} can be any of the
#' distributions in \code{timedist}. Similarly, \code{d____pt} and \code{p____pt} correspond
#' to the improper density and probability functions for promotion time models.
#'
#' The relationship between the linear predictor in \code{formula} and the time-to-event of
#' the non-cured elements is logarithmic as in the accelerated failure time models
#' (Lawless, 2003). However, for \code{cureformula}, if \code{ncausedist = "bernoulli"}
#' the relationship is similar to a logistic regression model
#' (\code{family = binomial(link = " logit ")} in \code{glm}) and if
#' \code{ncausedist = "poisson"} is similar to a Poisson model with logarithmic link
#' function (\code{family = poisson(link = "log")} in \code{glm}). See the references for
#' more details.
#'
#' @return A list of class "curereg" containing information about the fitted model.
#' Components of interest to users may include:
#'
#' \item{call }{
#' the matched call.
#' }
#' \item{coefficients }{
#' a named vector of coefficients obtained via Maximum Likelihood (see Details).
#' }
#' \item{std.error }{
#' a named vector of the estimated standard errors for the coefficients (see Details).
#' }
#' \item{vcov }{
#' A matrix of the estimated covariances between the coefficient estimates in the predictor
#' of the model.
#' }
#' \item{loglik }{
#' log-likelihood.
#' }
#' \item{AIC }{
#' AIC the (generalized) Akaike Information Criterion for the fitted model.
#' }
#'
#' @references Jackson, C.H. Christopher Jackson (2015). flexsurv: Flexible Parametric Survival and Multi-State Models.
#' R package version 0.7. https://CRAN.R-project.org/package=flexsurv.
#'
#' Maller, R. A., & Zhou, X. (1996). Survival analysis with long-term survivors. New York: Wiley.
#'
#' Lawless, J. F. (2011). Statistical models and methods for lifetime data (Vol. 362). John Wiley & Sons.
#'
#' Ortega, E. M., Cancho, V. G., & Paula, G. A. (2009). Generalized log-gamma regression models
#' with cure fraction. Lifetime Data Analysis, 15(1), 79-106.
#'
#' Peng, Y., Dear, K. B., & Denham, J. W. (1998). A generalized F mixture model for cure rate
#' estimation. Statistics in medicine, 17(8), 813-830.
#'
#' @author Rumenick Pereira da Silva <rumenickbf@hotmail.com>
#'
#' @seealso \code{\link{confint.curereg}} for confidence intervals for the coefficients, \code{\link{predict.curereg}} for predict cure rates
#' from fitted \code{curereg} model, \code{\link{plot.curereg}} and \code{\link{lines.curereg}} to
#' plot fitted survival, hazards and cumulative hazards from models fitted by \code{\link{curereg}}.
#'
#' @examples
#' ## fit Marshall-Olkin extended extreme value standart mixture model
#' data(e1684)
#' fitmo <- curereg(Surv(FAILTIME, FAILCENS) ~ TRT + SEX + AGE, cureformula = ~ TRT + SEX + AGE,
#'                 data = e1684, timedist = "moeev", ncausedist = "bernoulli")
#'
#' # Output of 'curereg' object
#' fitmo
#'
#' # Extract Model Coefficients
#' coef(fitmo)
#'
#' # Extract model coefficients:
#' # Terms: failure time distribution model
#' coef(fitmo, terms = "time")
#' # Terms: cure probability model
#' coef(fitmo, terms = "cure")
#'
#' # Object summaries
#' summary(fitmo)
#'
#' # Information criterion
#' AIC(fitmo) # Akaike information criterion
#' BIC(fitmo) # Bayesian information criterion
#'
#' # Extract Log-Likelihood
#' logLik(fitmo)
#' # Calculate Variance-Covariance Matrix for a Fitted Model Object
#' vcov(fitmo)
#'
#' @keywords models survival
#'
#' @importFrom  flexsurv flexsurvreg
#' @export
curereg <- function(formula, cureformula=~1, data, weights, timedist = "moeev",
                    ncausedist = "poisson", subset, na.action, inits, method = "SANN", ...) {

  call <- match.call()

  # Arguments for flexsurvreg function:
  argsfun <- list(...)
  if (is.null(cureformula) || missing(cureformula)) {
    dist <- timedist
    argsfun$anc <- NULL
    cureformula <- NULL
    dcure <- NULL
  }
  else {
    ncausedist <- match.arg(tolower(ncausedist), c("bernoulli", "poisson"))
    dcure <- ifelse(ncausedist=="bernoulli", "sm", "pt")
    dist <- match.arg(paste0(timedist, dcure), names(curereg.dists))
    argsfun$anc <- list(theta=cureformula)
    if (missing(inits)) {
      argsfun$inits <- NULL
      initp <- ifelse(missing(data),
                      1-mean(get(all.vars(formula)[2])),
                      ifelse(inherits(data, c("list", "data.frame")),
                             1-mean(data[[all.vars(formula)[2]]]),
                             stop("Must be data.frame or list")))
      curereg.dists[[dist]]$inits <- curereg.dists[[dist]]$inits(initp=initp)
    }
    else {
      argsfun$inits <- inits
      }
  }

  argsfun$control$maxit <- if (is.null(argsfun$control$maxit)) 10000 else argsfun$control$maxit
  argsfun$formula <- formula
  argsfun$data <- data
  argsfun$weights <- if(missing(weights)) NULL else weights
  argsfun$dist <- if (dist %in% names(curereg.dists)) curereg.dists[[dist]]  else dist
  argsfun$method <- method
  argsfun$subset <- if(missing(subset)) NULL else subset
  argsfun$na.action <- if(missing(na.action)) NULL else na.action

  # Fit model

  fit <- do.call(flexsurv::flexsurvreg, argsfun)

  # Output

  out <- list(call = call)
  out$formula <- formula
  out$cureformula <- cureformula

  if(!is.null(dcure)) out$dcure <- ifelse(dcure == "sm", "standart mixture model", "promotion time model") else out$dcure <- NULL
  out$X <- model.matrix(fit, fit$dlist$location)
  out$Z <- model.matrix(fit, "theta")

  if(is.null(names(fit$coef))) names.pars <- fit$dlist$par else names.pars <- names(fit$coef)

  for (i in seq_along(names.pars)) {
    if (any(names.pars[i] == c("sigma", "alpha", "shape", "P", "sdlog", "s1", "s2", "k"))) {
      names.pars[i] <- paste0("Log", "(", names.pars[i], ")")
    } else names.pars[i] <- names.pars[i]
  }

  if (is.null(argsfun$anc)) {
    if (dim(out$X)[2] == 1) {
      iloc <- grep(fit$dlist$location,names.pars)
      iother <- seq_along(names.pars)[-iloc]
      names.pars[iloc] <- "(Intercept)"
      names(fit$coef) <- names.pars
      ipar <- c(iloc, iother)
      out$coefficients <- list(coef.cure = NA, coef.time = fit$coef[ipar])
      out$std.error <- list(se.cure = NA, se.time = fit$res.t[ipar,"se"])
      names(out$std.error$se.time) <- names.pars[ipar]
      names(out$std.error$se.time)[1] <- paste(fit$dlist$location, colnames(out$X), sep=":")
      namesdim <- names(out$std.error$se.time)
    } else {
      iloc <-  which(names.pars %in% fit$dlist$location)
      ivar <- which(!(names(fit$coef) %in% fit$dlist$pars))
      iother <- seq_along(names.pars)[c(-iloc, -ivar)]
      ipar <- c(iloc, ivar, iother)
      names.pars[iloc] <- "(Intercept)"
      names(fit$coef) <- names.pars
      out$coefficients <- list(coef.cure = NA, coef.time = fit$coef[ipar])
      out$std.error <- list(se.cure = NA, se.time = fit$res.t[ipar,"se"])
      names(out$std.error$se.time) <- names.pars[ipar]
      names(out$std.error$se.time)[1:(length(ivar)+1)] <- paste(fit$dlist$location, colnames(out$X), sep=":")
      namesdim <- names(out$std.error$se.time)
    }
  } else {
    icure <- grep("theta", names.pars)
    iloc <- grep(fit$dlist$location, names.pars)
    itime <- which(names.pars %in% colnames(out$X))
    iother <- seq_along(names.pars)[-c(icure, iloc, itime)]
    ipartime <- c(iloc, itime, iother)
    ipar <- c(icure, ipartime)
    names(fit$coef) <- names.pars
    out$coefficients <- list(coef.cure = fit$coef[icure], coef.time = fit$coef[ipartime])
    out$std.error <- list(se.cure = fit$res.t[icure,"se"], se.time = fit$res.t[ipartime, "se"])
    names(out$coefficients$coef.time)[1:(length(itime)+1)] <- colnames(out$X)
    names(out$coefficients$coef.cure)[1:length(icure)] <- colnames(out$Z)
    names(out$std.error$se.time) <- names.pars[ipartime]
    names(out$std.error$se.time) [1:(length(itime)+1)] <- paste(fit$dlist$location, colnames(out$X), sep=":")
    names(out$std.error$se.cure) <- paste("theta", colnames(out$Z), sep=":")
    namesdim <- c(names(out$std.error$se.cure), names(out$std.error$se.time))
  }


  if (!is.matrix(fit$cov)) {
    out$vcov <- NA
  } else {
    if (!is.null(argsfun$fixedpars)) {
      ntpars <- length(ipar)
      out$vcov <- matrix(ncol = ntpars, nrow = ntpars)
      out$vcov[-argsfun$fixedpars, -argsfun$fixedpars] <- fit$cov
      aux <- numeric(length(argsfun$fixedpars))
      for(i in seq_along(argsfun$fixedpars)){ aux[i] <- which(argsfun$fixedpars[i] == ipar) }
      out$vcov <- out$vcov[ipar[-aux], ipar[-aux]]
      dimnames(out$vcov) <- list(namesdim[-aux], namesdim[-aux])
    } else {
      out$vcov <- fit$cov[ipar,ipar, drop = F]
      dimnames(out$vcov) <- list(namesdim, namesdim)
    }
  }

  out$npars <- fit$npars
  out$loglik <- fit$loglik
  out$AIC <- fit$AIC
  out$data <- model.frame(fit)
  out$opt <- fit$opt
  out$flexsurv <- fit
  rm(fit)
  out$ncausedist <- ncausedist
  out$timedist <- timedist
  class(out) <- "curereg"
  return(out)
}

#' @export
print.curereg  <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  distname <- switch(x$timedist,
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

  cat("Distribution:", paste(distname, ifelse(is.null(x$dcure), "", x$dcure), sep=" "), "\n")
  args <- list(...)
  if (is.null(args$digits)) args$digits <- 4
  if (is.null(args$scientific)) args$scientific <- FALSE
  cat("Coefficients:\nCure probability model:\n")
  if (is.null(x$cureformula)){
    cat("Not stated 'cureformula'!\n")
  } else{
    print(x$coefficients$coef.cure, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
    cat("\n")
  }

  cat("Failure time distribution model:\n")
  print(x$coefficients$coef.time, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
  n <- dim(x$data)[1]
  nevent <- sum(x$data[,1][,2])
  cat("\nn = ", n, ",", " Events: ", nevent, ",", " Censored: ", n-nevent, "\n",
      "Log-likelihood = ", x$loglik, "\nAIC = ", x$AIC, sep = "")
  cat("\n")
  invisible(x)
}

#' @export
summary.curereg  <- function(object, ...) {
  cat("Call:\n")
  dput(object$call)
  distname <- switch(object$timedist,
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

  cat("Distribution:", paste(distname, ifelse(is.null(object$dcure), "", object$dcure), sep=" "), "\n")
  args <- list(...)
  if (is.null(args$digits)) args$digits <- 4
  if (is.null(args$scientific)) args$scientific <- FALSE
  cat("\nCure probability model:\n")
  if (is.null(object$cureformula)){
    cat("Not stated 'cureformula'!\n")
  } else{
    zvalue <- object$coefficients$coef.cure/object$std.error$se.cure
    outcure <- cbind(object$coefficients$coef.cure, object$std.error$se.cure, zvalue, 1-pnorm(abs(zvalue)))
    colnames(outcure) <- c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)")
    print(outcure, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
    cat("\n")
  }

  cat("Failure time distribution model:\n")
  zvalue <- object$coefficients$coef.time/object$std.error$se.time
  outtime <- cbind(object$coefficients$coef.time, object$std.error$se.time, zvalue, 1-pnorm(abs(zvalue)))
  colnames(outtime) <- c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)")
  print(outtime, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
  n <- dim(object$data)[1]
  nevent <- sum(object$data[,1][,2])
  cat("\nn = ", n, ",", " Events: ", nevent, ",", " Censored: ", n-nevent, "\n",
      "Log-likelihood = ", object$loglik, "\nAIC = ", object$AIC, sep = "")
  cat("\n")
  invisible(object)
}

# terms = time or cure or missing (default)
#' @export
coef.curereg <- function(object, ...){
  if(missing(terms)) object$coefficients else object$coefficients[[paste0("coef.",terms)]]
}

# Predict cure fraction
#' @export
curefraction <- function(x, newData = NULL, unique = TRUE, ordered = TRUE, n = 6){
  if(!inherits(x, "curereg")) stop("Object must be results of curereg")
  curefun <- if(x$ncausedist == "poisson") function(x) exp(-exp(x)) else function(x) logit(x, inverse = TRUE)
  df <- if(is.null(newData)) x$Z else newData
  df <- if(unique) unique(df) else df
  cure <- curefun(tcrossprod(coef(x, terms = "cure"), df))[,]
  df <- if(ordered) cbind(df[,-1], curefraction = cure)[order(cure),] else cbind(df[,-1], curefraction = cure)
  cat("head:\n")
  print(head(df, n = n))
  cat("...\ntail:\n")
  print(tail(df, n = n))
  invisible(df)
}


# predict.curereg <- function(x, newX, newZ, ...) {
#
# }

#' @importFrom  flexsurv plot.flexsurvreg
#' @export
plot.curereg <- function(x, ...) {
  flexsurv::plot.flexsurvreg(x$flexsurv, ...)
}

#' @importFrom  flexsurv lines.flexsurvreg
#' @export
lines.curereg <- function(x, ...) {
  flexsurv::lines.flexsurvreg(x$flexsurv, ...)
}

#' @export
logLik.curereg <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$npars
  attr(val, "nobs") <-  nrow(object$data)
  class(val) <- "logLik"
  val
}

#' @export
vcov.curereg <- function(object, ...) {
  object$vcov
}

