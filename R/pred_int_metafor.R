#' Calculate prediction intervals using the 'metafor' pacakge
#' 
#' @title Prediction interval for a data.frame assuming effect size 'ROM' and 'metafor' package
#' @name pred_int_metafor
#' @description Prediction interval for a data.frame assuming effect size 'ROM'
#' @param x should be an object of class "data.frame"
#' @param method method used in computing the between-trial variance (default = "REML")
#' @param interval type of interval with default being "prediction"
#' @param var.names variable names for 'treatment', 'control' and 'id': c("trt","ctrl","ID") 
#' @param level coverage level with default 0.95
#' @return list with elements 'pdi' (prediction interval) and 'I2' (heterogeneity ratio)
#' @details see functions in metafor package 'escalc', 'rma' and 'predict.rma'
#' @export
#' @examples 
#' \dontrun{
#' ## Using soybean row spacing
#' data(soyrs)
#' pdi.I2 <- pred_int_metafor(soyrs, var.names = c("TRT1_Yld","TRT2_Yld","Trial_ID")) 
#' }
#'
#'

pred_int_metafor <- function(x, method = "REML", 
                        interval = c("prediction", "confidence"),
                        var.names = c("trt","ctrl","ID"),
                        level = 0.95){
  
  interval <- match.arg(interval)
  if(class(x) != "data.frame") stop("only for data frames")
  ## Make sure names align
  x2 <- x[,var.names]
  if(ncol(x2) != 3) stop("var.names might be wrong")
  ## Calculate preliminaries
  x.trt <- x[,var.names[1]]
  x.ctr <- x[,var.names[2]]
  trial <- x[,var.names[3]]
  x.m1i <- aggregate(x.trt ~ trial, data = x, FUN = mean)[,2]
  x.m2i <- aggregate(x.ctr ~ trial, data = x, FUN = mean)[,2]
  x.sd1i <- aggregate(x.trt ~ trial, data = x, FUN = sd)[,2]
  x.sd2i <- aggregate(x.ctr ~ trial, data = x, FUN = sd)[,2]
  x.n1i <- aggregate(x.trt ~ trial, data = x, FUN = length)[,2]
  x.n2i <- aggregate(x.ctr ~ trial, data = x, FUN = length)[,2]
  
  rr.escl <- escalc("ROM", m1i = x.m1i, m2i = x.m2i, 
                    sd1i = x.sd1i, sd2i = x.sd2i,
                    n1i = x.n1i, n2i = x.n2i)
  
  x.rma <- rma(rr.escl$yi, rr.escl$vi, data = rr.escl, method = method)
  
  x.pred <- unclass(predict(x.rma, level = level))
  if(interval == "prediction"){
    pdi <- c(x.pred$pred, x.pred$cr.lb, x.pred$cr.ub)
  }else{
    pdi <- c(x.pred$pred, x.pred$ci.lb, x.pred$ci.ub)
  }
  names(pdi) <- c("fit","lwr","upr")
  ans <- list(pdi = pdi, I2 = x.rma$I2)
  return(ans)
}