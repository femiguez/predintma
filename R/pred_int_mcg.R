#' Calculate prediction intervals using the 'MCMCglmm' pacakge and method 'ntrial'
#' 
#' @title Prediction interval for a data.frame using package 'MCMCglmm' and method 'ntrial'
#' @name pred_int_mcg_ntrial
#' @description Prediction interval for a data.frame using package 'MCMCglmm' and method 'ntrial'
#' @param formula an object of class formula with response ~ group
#' @param data a data frame which contains the response and the grouping 
#' @param level coverage level with default 0.95
#' @param control list which provides crontrol over the MCMC method \cr
#'                burnin = 5000 (default) \cr
#'                nitt = 50000 (default) \cr
#'                prior = list (see MCMCglmm) \cr
#' @return a prediction interval for a "new_trial"
#' @details see function 'MCMCglmm' in package 'MCMCglmm'. The method adds an extra row to 
#' the data.frame with missing response and then calculates the 'prediciton' for that new row
#' using the quantiles from the posterior distribution.
#' @export
#' @examples 
#' \dontrun{
#' ## Using soybean row spacing
#' data(soyrs)
#' pdi <- pred_int_mcg_ntrial(lrr ~ Trial_ID, data = soyrs)
#' pdi 
#' }
#'
#'

pred_int_mcg_ntrial <- function(formula, data, level = 0.95, 
                                control = list()){
  
  if(missing(data)) stop("data are missing")
  if(class(data) != "data.frame") stop("only for data frames")
  if(length(as.formula(formula)) != 3) stop("formula should be 'response' ~ 'trial'") 
  
  var.names <- all.vars(as.formula(formula))
  x <- subset(data, select = var.names)
  resp.name <- var.names[1]
  trial.name <- var.names[2]
  
  ## Set up alpha
  alpha <- 1 - level
  half.alpha <- alpha/2

  ## Create "new trial"
  ndat <- data.frame(x[1,])
  ndat[,trial.name] <- "A_new_trial"
  ndat[,resp.name] <- NA
  ## Created the data.frame with the "empty" trial
  dat <- rbind(x, ndat)
  
  mcmc.c <- mcmc.control()
  mcmc.c[names(control)] <- control
  prior1 <- mcmc.c$prior
  brn <- mcmc.c$burnin
  nitt <- mcmc.c$nitt
  
  trial.fm <- as.formula(paste0("~",trial.name))
  resp.fm <- as.formula(paste0(resp.name,"~1"))
  ## Run the chains
  mc1 <- mcparallel(MCMCglmm(fixed = resp.fm,
                   random = trial.fm, 
                   prior = prior1, pr = TRUE, 
                   nitt = nitt, burnin = brn,
                   data = dat, verbose = FALSE))
  
  mc2 <- mcparallel(MCMCglmm(fixed = resp.fm, 
                   random = trial.fm, 
                   prior = prior1, pr = TRUE, 
                   nitt = nitt, burnin = brn,
                   data = dat, verbose = FALSE))
  
  mcg1 <- mccollect(mc1)[[1]]
  mcg2 <- mccollect(mc2)[[1]]
  
  gd1 <- gelman.diag(mcmc.list(mcg1$Sol, mcg2$Sol))
  gd2 <- gelman.diag(mcmc.list(mcg1$VCV, mcg2$VCV))
  
  if(gd1$mpsrf > 1.1) warning("chain might have not converged, gd1 (fixed)")
  if(gd2$mpsrf > 1.1) warning("chain might have not converged, gd2 (random)")
  
  sdf <- as.data.frame(mcg2$Sol)
  ntn <- paste0(trial.name,".A_new_trial")
  new.t.beff <- sdf[,ntn]
  
  pd.mu <- as.data.frame(mcg2$Sol[,1])
  prdi.mu <- pd.mu + new.t.beff
  pdi <- quantile(prdi.mu[,1], probs = c(0.5, half.alpha, 1 - half.alpha))
  
  ans <- pdi
  
  return(ans)
}

#' @title Prediction interval for an object of class 'MCMCglmm' using method 'tdist'
#' @name pred_int_mcg_tdist
#' @description Prediction interval for an object of class 'MCMCglmm' using method 'tdist'
#' @param x should be an object of class 'MCMCglmm'
#' @param level coverage level with default 0.95
#' @return a prediction interval
#' @details see function 'MCMCglmm' in package 'MCMCglmm'. This method extracts the point estimates
#' for the mean, between-trial variance and within trial variance and plugs them in the 'tdist'
#' prediction formula 
#' @export

pred_int_mcg_tdist <- function(x, level = 0.95){
  ## This function extracts elements from a MCMCglmm object
  ## and calculates a prediction interval using the t-dist 
  ## approach
  if(class(x) != "MCMCglmm") stop("only for MCMCglmm objects")
  
  n.k <- length(x$Z@Dimnames[[2]])
  alph <- (1 - level)/2
  t.val <- qt(1 - alph, n.k-2) ## t.value
  mu <- summary(x)$solution[1]
  pvm <- summary(x$VCV)[[1]]
  tau.sqrd <- pvm[1] ## between studies variance
  if(nrow(summary(x$Sol)[[1]]) == 1){
    se.mu <- summary(x$Sol)[[1]][2] ## Intercept se
  }else{
    se.mu <- summary(x$Sol)[[1]][1,2] ## Intercept se
  }
  se.muh <- sqrt(tau.sqrd + se.mu^2) 
  pdi <- t.val * se.muh
  mu.lci <- mu - pdi
  mu.uci <- mu + pdi
  ans <- as.vector(c(mu, mu.lci, mu.uci))
  return(ans)
}

#' @name mcmc.control
#' @param burnin warmup or burnin iterations
#' @param nitt number of iterations for the MCMC method
#' @param prior prior specification as in the MCMCglmm package
#' @export
#' 
mcmc.control <- function(burnin = 5e3, nitt = 5e4, prior = NULL){
  
  if(missing(prior)){
    prior1 <- list(B = list(mu = 0, V = 10),
                   G = list(G1 = list(V = 1, nu = 0.002)),
                   R = list(V = 1, nu = 0.002))
  }
  ans <- list(burnin = burnin, nitt = nitt, prior = prior)
  ans
}