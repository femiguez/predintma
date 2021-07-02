
conf_int_stanreg <- function(x, level = 0.95, pmethod = c("epred")){
  
  if(!inherits(x, "stanreg")){
    stop("Object should be of class 'stanreg' ")
  }
  
  prbs1 <- (1 - level)/2
  prbs2 <- 1 - prbs1
  cnfi <- brms::posterior_summary(posterior_epred(x, re.form = NA), probs = c(prbs1, prbs2))
  ans <- colMeans(cnfi)
  return(ans[c(1,3,4)])
}

pred_int_stanreg <- function(x, level = 0.95, pmethod = c("epred", "predict","ntrial")){
  
  if(!inherits(x, "stanreg")){
    stop("Object should be of class 'stanreg' ")
  }
  
  pmethod <- match.arg(pmethod)
  
  prbs1 <- (1 - level)/2
  prbs2 <- 1 - prbs1
  ## This computes conditional modes (or 'BLUPs') and their distribution
  if(pmethod == "epred"){
    pep <- rstanarm::posterior_epred(x, re.form = NULL)  
  }
  if(pmethod == "predict"){
    pep <- rstanarm::posterior_predict(x, re.form = NULL)
  }
  if(pmethod == "ntrial"){
    ndat <- data.frame(Trial_ID = "new-trial")
    names(ndat) <- all.vars(x$formula)[2]
    pep <- rstanarm::posterior_epred(x, re.form = NULL, newdata = ndat)
  }
  
  pdi <- quantile(c(pep), probs = c(0.5, prbs1, prbs2))
  pdi
}