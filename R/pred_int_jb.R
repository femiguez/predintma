#' @title Prediction Interval for Random-effects Meta-analysis, which uses a combination of Jackknife and Bootstrap
#' @name pred_int_jb
#' @description Prediction for Random-effects Model using Subsampling
#' @param formula standard formula with response ~ trial
#' @param data should be an object of class "data.frame"
#' @param level coverage level with default 0.95
#' @param iter number of iterations for the bootstrap level 1
#' @param R number of iterations for the bootstrap level 2
#' @param var.names variable names to be passed to the 'data.frame' methods
#' @return a list with 'pdi' prediction intervals, 'cf' correction factor and 'coverage' empirical coverage
#' @export
#' @examples 
#' \dontrun{
#' ## Using soybean row spacing
#' data(soyrs)
#' pdi.rs <- pred_int_jb(formula = lrr ~ Trial_ID, data = soyrs) 
#' data(soysff)
#' pdi.ff <- pred_int_jb(formula = lrr ~ Trial_ID, data = soysff) 
#' }
#'

pred_int_jb <- function(formula, data, level = 0.95, iter = 10, R = 50, var.names=NULL){
  
  ## First step is I leave one study out
  
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  ## How many trials?
  trials <- unique(data[,var.names[2]])
  n.trials <- length(trials)
  x <- data[,var.names]
  names(x) <- c("y","trial")
  ## Need to redefine the individual level
  pdis <- data.frame(m = rep(NA, iter), lb = NA, ub = NA)
  cfs <- numeric(iter)
  curr.covs <- numeric(iter)
  alpha <- 1 - level
  half.alpha <- alpha/2
  cf <- 1
  incs <- numeric(iter)
  
  for(i in 1:iter){

    inc <- numeric(n.trials)
    for(j in 1:n.trials){
      ## Remove one of the trials
      xm1 <- subset(x,trial != trials[j])
      loo <- mean(subset(x, trial == trials[j])$y)
      ## Bootstrap means
      btm <- boot_tmeans(y ~ trial, data = xm1, R = R)
      ## Quantile prediction interval
      pdi.q <- quantile(btm$dat$ys, probs = c(0.5, half.alpha, 1 - half.alpha))
      ## Does it include the one that we left out?
      inc[j] <- includes(loo, pdi.q[2], pdi.q[3])
    }
    lb0 <- cf * abs(pdi.q[2] - pdi.q[1])
    lb <- pdi.q[1] - lb0
    ub0 <- cf * abs(pdi.q[3] - pdi.q[1])
    ub <- pdi.q[1] + ub0
    pdi <- c(pdi.q[1], lb, ub)
    pdis[i,1:3] <- as.vector(pdi)
    curr.cov <- mean(inc)
    curr.covs[i] <- curr.cov
    cf <- level / curr.cov
    cfs[i] <- cf
  }
 ans <- list(pdi=pdis, cf = cfs, coverage = curr.covs)
}