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

#' Calculate prediction interval based on the t-distribution. It uses 
#' 'pred_int_tdist_df'
#' 
#' 
#' @title Prediction interval for a 'data.frame' based on the t-distribution
#' @name pred_int_tdist_df
#' @description Simple wrapper which 'bootstraps' the trial means.
#' @param x should be an object of class "data.frame"
#' @param level coverage level with default 0.95
#' @param neval number of evaluations (default = 500)
#' @param interval type of interval, by default is "prediction"
#' @return a prediction interval
#' @seealso function 'pred_int_tdist_df'
#' @details the prediction is at the mean of 'x'. Missing values are not allowed.
#' @export
#' @examples 
#' \dontrun{
#' ## Prediciton interval for soybean data
#' data(soyrs)
#' pdi <- pred_int_tdist_df(soyrs, var.names = c("lrr","Trial_ID")) 
#' }
#'
pred_int_tdist_df <- function(x, level = 0.95, neval = 500, 
                              interval = c("prediction","confidence"),
                              var.names = c("y","trial")){
  
  interval <- match.arg(interval)
  ## How many trials?
  trials <- unique(x[,var.names[2]])
  n.trials <- length(trials)
  ## Need to redefine the individual level
  alpha <- 1 - level
  i.level <- 1 - alpha/neval
  ## set up dataset
  x2 <- x[,var.names]
  names(x2) <- c("y","trial")
  ## do this 'neval' number of times
  pdi.c <- matrix(NA, ncol = 3, nrow = neval)
  for(i in 1:neval){
    sub.sample <- numeric(n.trials)
    k <- 1
    for(j in trials){
      xs <- subset(x2, trial == j)[,"y"]
      sub.sample[k] <- sample(xs, size = 1)
      k <- k + 1
    }
    pdi.c[i,] <- pred_int_tdist(sub.sample, level = level, interval = interval)
  }
  ans <- colMeans(pdi.c)
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
#' @param point.pred point prediction either 'mean' or 'median'
#' @param padding how much padding to add to the linear evaluation (default 0.0)
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
#' 
#' ## Following example from the link above
#' set.seed(12345)
#' n <- 1000
#' U <- rnorm(n)
#' hist(U)
#' u <- 2 ## New proposed value
#' Ua <- c(U,u) ## Augmented set
#' hist(Ua)
#' abline(v=u)
#' print(round(mean(Ua > u), 3))
#' nEval <- 200
#' u.candidate <- seq(from=min(U), to=max(U), length = nEval)
#' Cbounds <- c(u.candidate[match(5, apply(as.matrix(u.candidate),
#'                                         1,
#'                                        function(x){floor(2.5+100*(mean(c(U, x) < x)))}))],
#'             rev(u.candidate)[match(95, rev(apply(as.matrix(u.candidate),
#'                                                  1,
#'                                                  function(x){ceiling(-2.5+100*(mean(c(U, x) < x)))})))])
#' print(Cbounds)
#' hist(U)
#' abline(v = Cbounds)
#' ## Alternatively
#' pdi.cb <- pred_int_conformal(U)
#' 
#' ## From Shafer and Vovk (pg. 5)
#' ## Numners generated by Emmanuel Czuber
#' xn <- c(17,20,10,17,12,15,19,22,17,19,14,22,18,17,13,12,18,15,17)
#' pdi.xn <- pred_int_conformal(xn, neval = 500)
#' }
#'
#'
#'
pred_int_conformal <- function(x, neval = 200, level = 0.95,
                               point.pred = c("mean","median"),
                               padding = 0.0,
                               method = c("quantile","deviation")){
  
  if(class(x) != "numeric") stop("class should be numeric")
  point.pred <- match.arg(point.pred)
  method <- match.arg(method)
  
  ## This is our linear evaluation
  x10 <- max(x) * padding
  xx <- seq(from = min(x) - x10, to = max(x) + x10, length.out = neval)
  
  if(method == "quantile"){
    ## Lower percent
    alpha <- 1 - level
    lp <- (alpha / 2) * 100
  
    for(i in 1:neval){
      ## Augmented sample
      x.a <- c(x, xx[i])
      ## The '<' operation is logical and it prints
      ## TRUE (1) or FALSE (0)
      ## Calculate lower bound
      lb.x.a <- lp + 100 * mean(x.a < xx[i])
      if(lb.x.a >= (alpha * 100)){
        lower.bound <- xx[i]
        break
      }
    }
  
    for(i in neval:1){
      ## Augmented sample
      x.a <- c(x, xx[i])
      ub.x.a <- -lp + 100 * mean(x.a < xx[i])
      if(ub.x.a <= ((1 - alpha) * 100)){
        upper.bound <- xx[i]
        break
      }
    }
  }
  
  if(method == "deviation"){
    ## Calculate the sum of observations
    sum.x <- sum(x)
    length.p1 <- length(x) + 1
    cset <- numeric(neval)
    k <- 1
    for(i in 1:neval){
      tmp1 <- abs(sum.x - length(x) * xx[i])
      tmp2 <- abs(sum.x + xx[i] - (length.p1 * max(x)))
      tmp3 <- abs(sum.x + xx[i] - (length.p1 * min(x)))
      if(tmp1 <= max(c(tmp2,tmp3))){
        cset[k] <- xx[i]
        k <- k + 1
      }
    }
    cset <- cset[1:(k-1)]
    lower.bound <- min(cset)
    upper.bound <- max(cset)
  }
  
  if(point.pred == "median"){
    pprd <- median(x)
  }else{
    pprd <- mean(x)
  }
  ans <- c(pprd, lower.bound, upper.bound)
  return(ans)
}