#' Calculate nuisance parameter (degrees of freedom)
#' 
#' 
#' @title degrees of freedom calculation using 'emmeans'
#' @name npv
#' @description Calculates the nuisance parameter (degrees of freedom) using different methods
#' @param x should be an object of class "lmerMod"
#' @param degfr degrees of freedom method ("default","zdist","kr") 
#' @param level coverage level with default 0.95
#' @return a single value 'numeric' to be used as input to a t-distribution
#' @export
#' @examples 
#' \dontrun{
#' ## Include a meaningful example 
#' }
#'
#'

npv <- function(x, degfr = c("default","zdist","kr"), level = 0.95){
  
  degfr <- match.arg(degfr)
  if(class(x) != "lmerMod") stop("class should be 'lmerMod'")
  ## Calculate level for all of them 
  alph <- 1 - level
  qnt <- 1 - alph/2
  ## np stands for nuiscance parameter
  if(degfr == "default"){
    ## Number of trials
    nbrt <- getME(x, "q") ## or x@Gp[2]
    ## I substract 2 because this is in the
    ## original formula by Higgins, but it
    ## could be argued that I should subtract 
    ## 3 instead
    defr <- nbrt - 2
    npv <- qt(qnt, defr) 
  }
  
  if(degfr == "zdist"){
    ## t-dist with "Inf" degfr
    npv <- qnorm(qnt) 
  }
  
  if(degfr == "kr"){
    ## Kenward-Roger option
    degfr <- summary(emmeans(x, specs = ~1,
                             lmer.df = "kenward-roger"))$df
    npv <- qt(qnt, degfr) 
  }
  
##  if(degfr == "sat"){
##    stop("dont't want to use this method at the moment")
##    require(emmeans)
##    require(lmerTest)
##    ## Satterthwaite option
##    degfr <- summary(emmeans(x, specs = ~1, 
##                             lmer.df = "satterthwaite"))$df
##    npv <- qt(qnt, degfr) 
##  }
  
  npv
}

#' Calculate prediction interval based on the t-distribution. It uses 
#' 'lm' and 'predict'
#' 
#' 
#' @title Prediction interval for a 'numeric' based on the t-distribution
#' @name pred_int_tdist
#' @description Simple wrapper for 'lm' and 'predict' which returns a prediction interval. 
#' This method assumes iid normally distributed observations.
#' @param x should be an object of class "numeric"
#' @param level coverage level with default 0.95
#' @param interval type of interval, by default is "prediction"
#' @return a prediction interval
#' @seealso function 'lm' and 'predict.lm'
#' @details the prediction is at the mean of 'x'. Missing values are not allowed.
#' @export
#' @examples 
#' \dontrun{
#' ## Prediciton interval for simulated data
#' x <- rnorm(100)
#' pdi <- pred_int_tdist(x) 
#' }
#'
#'
#'
pred_int_tdist <- function(x, level = 0.95, interval = c("prediction","confidence")){
  
  interval <- match.arg(interval)
  if(class(x) != "numeric") stop("class should be numeric")
  
  ans <- predict(lm(x ~ 1), newdata=data.frame(x = mean(x)), 
                   interval = interval, level = level)
  return(ans)
}


#' Calculate prediction interval based on conformal inference. 
#' 
#' 
#' @title Prediction interval for a 'numeric' based on conformal inference
#' @name pred_int_conformal
#' @description Conformal inference prediction interval for a numeric
#' @param x should be an object of class "numeric"
#' @param neval number of evaluations to perform
#' @param level coverage level with default 0.95
#' @param point.pred point prediction either 'median' or 'mean'
#' @return a prediction interval
#' @details I wrote this function after reading this tutorial: https://cdsamii.github.io/cds-demos/conformal/conformal-tutorial.html
#' @export
#' @examples 
#' \dontrun{
#' ## Prediciton interval for simulated data
#' x <- rnorm(100)
#' pdi.t <- pred_int_tdist(x) 
#' pdi.c <- pred_int_conformal(x)
#' ## Not really a prediction interval, but related
#' pdi.q <- quantile(x, probs = c(0.5, 0.025, 0.975))
#' }
#'
#'
#'
pred_int_conformal <- function(x, neval = 200, level = 0.95,
                               point.pred = c("median","mean")){
  
  if(class(x) != "numeric") stop("class should be numeric")
  point.pred <- match.arg(point.pred)
  ## This is our linear evaluation
  xx <- seq(from = min(x), to = max(x), length.out = neval)
  
  ## Lower percent
  alpha <- 1 - level
  lp <- (alpha / 2) * 100
  
  for(i in 1:neval){
    ## Augmented sample
    x.a <- c(x, xx[i])
    ## The '<' operation is logical and it prints
    ## TRUE (1) or FALSE (0)
    ## Calculate lower bound
    lb.x.a <- floor(lp + 100 * mean(x.a < xx[i]))
    if(lb.x.a >= alpha * 100){
      lower.bound <- xx[i]
      break
    }
  }
  
  for(i in neval:1){
    ## Augmented sample
    x.a <- c(x, xx[i])
    ub.x.a <- ceiling(-lp + 100 * mean(x.a < xx[i]))
    if(ub.x.a <= (1 - alpha) * 100){
      upper.bound <- xx[i]
      break
    }
  }
  
  if(point.pred == "median"){
    pprd <- median(x)
  }else{
    pprd <- mean(x)
  }
  ans <- c(pprd, lower.bound, upper.bound)
  return(ans)
}