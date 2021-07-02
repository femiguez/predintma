#' 
#' 
#' @title Prediciton interval for Inferential Models
#' @name pred_int_im
#' @param x object of class \sQuote{merMod}
#' @param level probability level for interval coverage
#' @return a prediction interval
#' @export
#' @examples 
#' \donttest{
#' require(lme4)
#' data(soyrs)
#' fit <- lmer(lrr ~ 1 + (1|Trial_ID), data = soyrs)
#' pdi <- pred_int(fit)
#' pdi 
#' ## New formula
#' pdi2 <- pred_int_im(fit)
#' 
#' }

pred_int_im <- function(x, level = 0.95, dist = c("z", "t")){
  
  if(!inherits(x, "merMod"))
    stop("This function is for objects which inherit class 'merMod'")
  
  dist <- match.arg(dist)
  
  SE1 <- sigma(x)^2 * (1/nobs(x) + 1/ngrps(x))
  ## Number of samples per group
  n.k <- as.vector(table(x@flist[[1]]))^2
  SE2 <- VarCorr(x)[[1]][1] * (1 + 1/nobs(x)^2 * sum(n.k))
  SE <- sqrt(SE1 + SE1)
  
  if(dist == "z"){
    Zlwr <- stats::qnorm((1 - level)/2)  
  }else{
    Zlwr <- -npv(x, degfr = "default", level = level)
  }
  
  lb <- fixef(x) + Zlwr * SE
  ub <- fixef(x) - Zlwr * SE
  
  ans <- c(fixef(x), lb, ub)
  names(ans) <- c("fit", "lwr", "upr")
  
  return(ans)
}