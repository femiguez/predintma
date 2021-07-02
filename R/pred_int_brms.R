

conf_int_brmsfit <- function(x, level = 0.95){
  
  if(!inherits(x, "brmsfit")){
    stop("Object should be of class 'brmsfit' ")
  }
  
  prbs1 <- (1 - level)/2
  prbs2 <- 1 - prbs1
  cnfi <- posterior_summary(posterior_epred(x, re_formula = NA), probs = c(prbs1, prbs2))
  ans <- colMeans(cnfi)
  return(ans[c(1,3,4)])
}

pred_int_brmsfit <- function(x, level = 0.95, pmethod = c("epred", "predict")){
  
  if(!inherits(x, "brmsfit")){
    stop("Object should be of class 'brmsfit' ")
  }
  pmethod <- match.arg(pmethod)
  
  prbs1 <- (1 - level)/2
  prbs2 <- 1 - prbs1
  ## This computes conditional modes (or 'BLUPs') and their distribution
  if(pmethod == "epred"){
    pep <- posterior_epred(x, re_formula = NULL)  
  }else{
    pep <- posterior_predict(x, re_formula = NULL)
  }
  
  pdi <- quantile(c(pep), probs = c(0.5, prbs1, prbs2))
  pdi
}