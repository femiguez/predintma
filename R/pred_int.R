#' This function will calculate a prediction interval in the context of
#' meta-analysis
#' As defined in Higgins et al. (2009) it considers the parameteric uncertainty
#' and the between study uncertainty, but there is more to this...
#' Implemented cases:
#'
#' object of class 'numeric' \cr
#' object of class 'lmerMod' \cr
#' object of class 'MCMCglmm' \cr
#' object of class 'data.frame' \cr
#'
#' For object of class numeric method t-dist is based on the 'lm' function and predict
#' There is also a method based on conformal prediction. See 'pred_int_conformal'
#'
#' For object of class 'lmerMod' there are three methods: \cr
#' - tdist \cr
#' - tdist2 \cr
#' - boot (uses lme4::bootMer) \cr
#'
#' For object of class 'MCMCglmm' there are four methods: \cr
#' - tdist \cr
#' - mcmc \cr
#' - simulate \cr
#' - predict \cr
#' 
#' There is also the 'ntrial' method, which needs to
#' be applied to a data.frame. It uses 'MCMCglmm'
#' 
#' Another method 'pred_int_metafor' uses the metafor pacakge
#' And it also requires a data.frame
#' 
#' the degrees of freedom argument (degfr) allows
#' for different methods other than n.k-2
#' Calculate prediction intervals for a variety of objects
#' 
#' @title Prediction interval for different objects
#' @name pred_int
#' @description Prediction interval for a variety of objects in the context of random-effects meta-analysis
#' @param x should be an object of class: "numeric", "data.frame", "lmerMod" or "MCMCglmm"
#' @param interval type of interval with default 'prediction'
#' @param method either 'tdist', 'tdist2', 'boot', 'mcmc', 'simulate', 'predict',
#'               'metafor', 'ntrial' or 'conformal'.
#' @param m.method method for between-trial variance estimator (used in the 'metafor' package)
#' @param c.method method for conformal prediction. See ?pred_int_conformal
#' @param degfr degrees of freedom method default (n.k-2), zdist ("Inf") or "kr" (Kenward-Roger).
#'              see package 'emmeans'
#' @param level coverage level with default 0.95
#' @param nsim number of simulations for the 'boot' method.
#' @param var.names variable names to be passed to the 'data.frame' methods
#' @param formula formula interface for 'data.frame' methods
#' @param ... arguments to be passed to a few of the functions
#' @return a prediction interval for a "new_trial"
#' @details The main reference to understand the rationale behind this function is \cr
#' Higgins et al. (2009) A re-evaluation of random-effect meta-analysis
#' @export
#' @examples 
#' \dontrun{
#' ## Using soybean row spacing
#' data(soyrs)
#' fit <- lmer(lrr ~ 1 + (1|Trial_ID), data = soyrs)
#' pdi <- pred_int(fit)
#' pdi 
#' }
#'
#'

pred_int <- function(x,  
                     interval = c("prediction","confidence"),
                     method = c("tdist","tdist2","boot","mcmc","simulate","predict",
                                "metafor","ntrial","conformal"),
                     m.method = c("REML","DL", "HE", "SJ", "ML", "REML", "EB", "HS", "GENQ"),
                     c.method = c("quantile","deviation","jackknife"),
                     degfr = c("default","zdist","kr"),
                     level = 0.95,
                     nsim = 500,
                     var.names = NULL,
                     formula = NULL,...){
  
  interval <- match.arg(interval)
  method <- match.arg(method)
  m.method <- match.arg(m.method)
  degfr <- match.arg(degfr)
  
  if(!inherits(x, "numeric") && 
     !inherits(x, "lmerMod") &&
     !inherits(x, "MCMCglmm") &&
     !inherits(x, "data.frame")) stop("object not supported")
  
  if(inherits(x, "lmerMod") && interval == "confidence"){ 
    stop("not implemented, use 'confint' instead")
  }
  
  if(inherits(x, "numeric") && method == "tdist"){
    ans <- predict(lm(x ~ 1), newdata = data.frame(x = mean(x)), interval = interval, level = level)
  }
  
  if(inherits(x, "numeric") && method == "conformal"){
    ans <- pred_int_conformal(x, level = level, method = m.method, ...)
  }
  
  if(inherits(x, "data.frame") && method == "metafor"){
    ans <- pred_int_metafor(x, method = m.method, 
                            var.names = var.names,
                            level = level,
                            interval = interval)$pdi
  }
  
  if(inherits(x, "data.frame") && method == "ntrial"){
    ans <- pred_int_mcg_ntrial(formula = formula, data = x, level = level)
  }
  
  if(inherits(x, "lmerMod") && method == "tdist"){
    ## Extract needed components
    n.k <- x@Gp[2] ## Number of trials
    ## New method for calculating degrees of freedom
    ## npv = nuiscance parameter value
    t.val <- npv(x, degfr = degfr, level = level)
    ## old method 
    ## alph <- (1 - level)/2
    ## t.val <- qt(1 - alph, n.k-2) ## t.value
    mu <- fixef(x)
    tau.sqrd <- VarCorr(x)[[1]][1] ## between studies variance
    se.mu <- summary(x)$coefficients[2] ## Intercept standard error
    ## The following line combines parametric uncertainty (se.mu)
    ## With sampling uncertainty (tau.sqrd)
    se.muh <- sqrt(tau.sqrd + se.mu^2) 
    pdi <- t.val * se.muh
    mu.lci <- mu - pdi
    mu.uci <- mu + pdi
    ans <- as.vector(c(mu, mu.lci, mu.uci))
  }
  
  if(inherits(x, "lmerMod") && method == "tdist2"){
    ## Extract needed components
    n.k <- x@Gp[2] ## Number of trials
    ## Coming up with weights
    ## mean within trial sample size
    ## mni <- mean(colSums(getME(x, "Z")))
    ## w1 <- n.k/(n.k + mni)
    ## w2 <- mni/(n.k + mni)
    if(n.k < 3) stop("number of studies too small")
    ## New degrees of freedom adjustment
    ## with options
    t.val <- npv(x, degfr = degfr, level = level)
    mu <- fixef(x)
    tau.sqrd <- VarCorr(x)[[1]][1] ## between studies variance
    se.mu <- summary(x)$coefficients[2] ## Intercept standard error
    ## Miguez correction
    sigma.sqrd <- sigma(x)^2
    icc <- tau.sqrd / (tau.sqrd + sigma.sqrd)
    sigma.sqrd.c <- (1 - icc) * sigma.sqrd ## Correction
    ## The following line combines parametric uncertainty (se.mu)
    ## With sampling uncertainty (tau.sqrd)
    ## unweighted
    se.muh <- sqrt(sigma.sqrd.c + icc*tau.sqrd + se.mu^2) 
    ## weighted
    ## se.muh <- sqrt(w2*sigma.sqrd + w1*icc*tau.sqrd + se.mu^2) 
    pdi <- t.val * se.muh
    mu.lci <- mu - pdi
    mu.uci <- mu + pdi
    ans <- as.vector(c(mu, mu.lci, mu.uci))
  }
  
  if(inherits(x, "lmerMod") && method == "boot"){
    tmp <- suppressWarnings(bootMer(x, 
                                    nsim = nsim, 
                                    pred_int, 
                                    use.u = TRUE)$t)
    ans <- colMeans(tmp)
  }
 
  if(inherits(x, "MCMCglmm") && method == "tdist"){
    ans <- pred_int_mcg_tdist(x, level = level)
  }
  
  if(inherits(x, "MCMCglmm") && method == "mcmc"){
    mu.chain <- as.vector(x$Sol[,1])
    tau.chain <- sqrt(as.vector(x$VCV[,1]))
    err.chain <- rnorm(length(mu.chain), 0, tau.chain)
    prd.chain <- mu.chain + err.chain
    mu.pi <- c(median(prd.chain), quantile(prd.chain, probs = c(0.025,0.975)))
    ans <- mu.pi
  }
  
  if(inherits(x, "MCMCglmm") && method == "simulate"){
    ## Extract number of effective samples
    nrs <- nrow(as.matrix(x$Sol))
    ndat <- data.frame(lrr = rep(0,nrs), 
                       Trial_ID = "A_new_trial")
    prd.b <- scale(simulate(x, newdata = ndat,
                           marginal = NULL), scale = FALSE)
    pd.mu <- as.matrix(x$Sol[,1])
    prdi.c.mu <- pd.mu + prd.b
    ans <- quantile(prdi.c.mu, probs = c(0.5, 0.025,0.975))
  }
  
  if(inherits(x, "MCMCglmm") && method == "predict"){
    ## I will assume the pr = TRUE when the model was ran
    ## It is assumed that there is a factor called 'Trial_ID'
    if(as.character(x$Random$formula)[2] != "Trial_ID")
      stop("Trial_ID not found")
    tns <- x$Z@Dimnames[[2]]
    trnms <- sapply(tns, FUN = function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])
    rsp.nm <- as.character(x$Fixed$formula[2])
    ## The next line is the median for the intercept
    rsp.m <- as.data.frame(summary(x$Sol)[[2]][1,3])
    names(rsp.m) <- rsp.nm
    ndat <- data.frame(rsp.m, 
                       Trial_ID = as.vector(trnms))
    ans <- colMeans(predict(x, 
                            newdata = ndat, 
                            marginal = NULL, 
                            interval = interval))
  }
  
  setNames(ans, c("fit","lwr","upr"))
  return(ans)
}




