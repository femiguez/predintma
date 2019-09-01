
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
  if(class(x) != "numeric" & class(x) != "integer") stop("class should be numeric or integer")
  
  ans <- predict(lm(x ~ 1), newdata=data.frame(x = mean(x)), 
                 interval = interval, level = level)
  return(ans)
}

#' Calculate prediction interval based on the t-distribution. It uses 
#' 'pred_int_tdist_df'
#' 
#' @title Prediction interval for a 'data.frame' based on the t-distribution
#' @name pred_int_tdist_df
#' @description Simple wrapper which 'bootstraps' the trial means.
#' @param x should be an object of class "data.frame"
#' @param level coverage level with default 0.95
#' @param neval number of evaluations (default = 500)
#' @param interval type of interval, by default is "prediction"
#' @param var.names names of variables for response and 'trial'
#' @param formula alternative formula interface
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
                              var.names = c("y","trial"),
                              formula = NULL){
  
  interval <- match.arg(interval)
  
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  ## How many trials?
  trials <- unique(x[,var.names[2]])
  n.trials <- length(trials)
  ## Need to redefine the individual level
  alpha <- 1 - level
  i.level <- 1 - alpha/neval
  ## set up dataset
  x2 <- x[,var.names]
  if(dim(x2)[2] != 2) stop("var.names might be wrong")
  ## do this 'neval' number of times
  pdi.c <- matrix(NA, ncol = 3, nrow = neval)
  for(i in 1:neval){
    sub.sample <- numeric(n.trials)
    k <- 1
    for(j in trials){
      xs <- x2[x2[,2] == j, var.names[1]]
      sub.sample[k] <- sample(xs, size = 1)
      k <- k + 1
    }
    pdi.c[i,] <- pred_int_tdist(sub.sample, level = level, interval = interval)
  }
  ans <- colMeans(pdi.c)
  names(ans) <- c("fit","lwr","upr")
  return(ans)
}
