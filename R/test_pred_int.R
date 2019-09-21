#' Test the coverage for different prediction methods
#'
#' @title test_pred_int
#' @param formula standard R formula, 'response ~ trial'
#' @param data a data frame
#' @param method either 'tdist', 'ntrial', 'jags' or 'boot'
#' @param level prediction level
#' @param iter number of iterations
#' @param ... additional arguments to be passed to the functions
#' @export
#' @examples 
#' \dontrun{
#'  data(soyrs)
#'  t.tdist <- test_pred_int(lrr ~ Trial_ID, data = soyrs)
#' }
#'
#'
#'
#'

test_pred_int <- function(formula, data, 
                          method = c("tdist","ntrial","jags","boot"), 
                          level=0.95, iter = 1, quiet = TRUE, ...){
  
  ## First: extract variables
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  
  ## Method match
  method <- match.arg(method)
  
  dat <- data[,var.names]
  names(dat) <- c("y","trial")
  
  trials <- unique(dat[,"trial"])
  n.trials <- length(trials)
  
  pdis <- data.frame(m = rep(NA, iter), lb = NA, ub = NA)
  coverage <- numeric(iter)
  incs <- numeric(iter)
  
  start <- Sys.time()
  
  for(i in 1:iter){
    inc <- numeric(n.trials)
    for(j in 1:n.trials){
      ## Remove one of the trials
      xm1 <- droplevels(subset(dat, trial != trials[j]))
      loo <- mean(subset(dat, trial == trials[j])$y)
      ## Implement the different methods
      if(method == "tdist"){
        fit <- lmer(y ~ 1 + (1 | trial), data = xm1)
        pdi <- pred_int(fit, level = level)
      }
      if(method == "ntrial"){
        pdi <- pred_int_mcg_ntrial(formula = y ~ trial, data = xm1, level = level)
      }
      if(method == "jags"){
        pdi <- pred_int_jags(y ~ trial, data = xm1, level = level, quiet = TRUE, return = "pdi")
      }
      if(method == "boot"){
        pdi <- pred_int_boot(y ~ trial, data = xm1, level = level)
      }
      ## Does it include the one that we left out?
      inc[j] <- includes(loo, pdi[2], pdi[3])
    }
    pdis[i,1:3] <- as.vector(pdi)
    coverage[i] <- mean(inc)
    if(!quiet){
      if(i %% (n.trials%/%10) == 0){
        prog.time <- Sys.time() - start
        cat("progress ",i/n.trials*100,"% time: ",prog.time,"\n")
      }
    }
  }
  ans <- list(pdi=pdis, coverage = coverage)
  ans
}