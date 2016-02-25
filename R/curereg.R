############################################################################
#
# toInstall <- c("ggplot2", "reshape2", "RColorBrewer")
# doInstall <- toInstall %in% rownames(installed.packages())
# doInstall
# if(doInstall){install.packages(toInstall, repos = "http://cran.us.r-project.org")}
# lapply(toInstall, library, character.only = TRUE)

# Just-in-time compilation
compiler::enableJIT(3)

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

# list of avaliable distribuition:
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

#  Função para ajuste de modelos com fração de cura:

curereg <- function(formula, cureformula=~1, data, weights = NULL, timedist = "moeev",
                    ncausedist = "poisson", subset = NULL, na.action = "na.omit", inits = NULL,
                    method = "SANN", ...) {
   call <- match.call()
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
     initp <- ifelse(missing(data),
                     1-mean(get(all.vars(formula)[2])),
                     ifelse(inherits(data, c("list", "data.frame")),
                            1-mean(data[[all.vars(formula)[2]]]),
                            stop("Must be data.frame or list")))
     curereg.dists[[dist]]$inits <- curereg.dists[[dist]]$inits(initp=initp)
     argsfun$anc <- list(theta=cureformula)
     }

  # Argumentos a serem passados para função flexsurvreg
  argsfun$control$maxit <- if (is.null(argsfun$control$maxit)) 10000 else argsfun$control$maxit
  argsfun$formula <- formula
  argsfun$data <- data
  argsfun$weights <- weights
  argsfun$inits <- inits
  argsfun$dist <- if (dist %in% names(curereg.dists)) curereg.dists[[dist]]  else dist
  argsfun$method <- method
  argsfun$subset <- subset
  argsfun$na.action <- na.action
  # Ajustando modelo
  fit <- do.call(flexsurv::flexsurvreg, argsfun)
  # Construindo saída
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
         icure <- grep("theta",names.pars)
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
#   } else {
#     if(!is.null(argsfun$fixedpars)){ Se fixedpars = TRUE
#       aux <- numeric(length(argsfun$fixedpars))
#       for(i in seq_along(argsfun$fixedpars)){ aux[i] <- which(argsfun$fixedpars[i] == ipar) }
#       out$vcov <- fit$cov[ipar[-aux],ipar[-aux], drop=F]
#       dimnames(out$vcov) <- list(namesdim[-aux], namesdim[-aux])
     } else{
      out$vcov <- fit$cov[ipar,ipar, drop=F]
      dimnames(out$vcov) <- list(namesdim, namesdim)
    }
#  }

  out$npars <- fit$npars
  out$loglik <- fit$loglik
  out$AIC <- fit$AIC
  out$data <- model.frame(fit)
  out$opt <- fit$opt
  rm(fit)
  out$ncausedist <- ncausedist
  out$timedist <- timedist
  class(out) <- "curereg"
  return(out)
}

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

summary.curereg  <- function(x, ...) {
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
  cat("\nCure probability model:\n")
  if (is.null(x$cureformula)){
    cat("Not stated 'cureformula'!\n")
  } else{
    zvalue <- x$coefficients$coef.cure/x$std.error$se.cure
    outcure <- cbind(x$coefficients$coef.cure, x$std.error$se.cure, zvalue, 1-pnorm(abs(zvalue)))
    colnames(outcure) <- c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)")
    print(outcure, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
    cat("\n")
  }

  cat("Failure time distribution model:\n")
  zvalue <- x$coefficients$coef.time/x$std.error$se.time
  outtime <- cbind(x$coefficients$coef.time, x$std.error$se.time, zvalue, 1-pnorm(abs(zvalue)))
  colnames(outtime) <- c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)")
  print(outtime, print.gap=2, digits = args$digits, scientific = args$scientific, quote=FALSE, na.print="")
  n <- dim(x$data)[1]
  nevent <- sum(x$data[,1][,2])
  cat("\nn = ", n, ",", " Events: ", nevent, ",", " Censored: ", n-nevent, "\n",
      "Log-likelihood = ", x$loglik, "\nAIC = ", x$AIC, sep = "")
  cat("\n")
  invisible(x)
}

# terms = time or cure or NULL (default)
coef.curereg <- function(x, terms=NULL){
  if(is.null(terms)) x$coefficients else x$coefficients[[paste0("coef.",terms)]]
}

# Usar flexsurv::model.frame.flexsurvreg
predict.curereg <- function(x, newX, newZ, ...) {

}

# Predizer a fração de cudaros

curefraction <- function(x, newData = NULL, unique = TRUE, ordered = TRUE, n = 6){
  curefun <- if(x$ncausedist == "poisson") function(x) exp(-exp(x)) else function(x) logit(x, inverse = T)
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


plot.curereg <- function(x, ...) {

}


lines.curereg <- function(x, ...) {

}

logLik.curereg <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$npars
  attr(val, "nobs") <-  nrow(object$data)
  class(val) <- "logLik"
  val
}

vcov.curereg <- function(object, ...) {
  object$vcov
}

