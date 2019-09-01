#' @title Conformal Prediction Interval for Random-effects Meta-analysis
#' @name pred_int_conformal_df
#' @description Conformal Prediction for Random-effects Model using Subsampling
#' @param x should be an object of class "data.frame"
#' @param method "subsampling"
#' @param level coverage level with default 0.95
#' @param neval.1 number of evaluations at the 'trial' level.
#' @param neval.2 number of evaluations at the 'rep' level.
#' @param s.size within-group sample size, when the number of groups is less than 20, it should be higher than 1. 
#' @param replace replacement option in the within-group sampling
#' @param point.pred 'point prediction', either 'mean' or 'median'
#' @param var.names variable names to be passed to the 'data.frame' methods
#' @param formula alternative formula interface, if this is used, 'var.names' is ignored
#' @return a prediction interval for a "new_trial"
#' @details TODO
#' @export
#' @examples 
#' \dontrun{
#' ## Using soybean row spacing
#' data(soyrs)
#' pdi.conf <- pred_int_conformal_df(formula = lrr ~ Trial_ID, x = soyrs) 
#' }
#'
#'
#'

pred_int_conformal_df <- function(x, method = "subsampling", level = 0.95, 
                                  neval.1 = 500, neval.2 = 500, 
                                  s.size = 1, replace = FALSE,
                                  point.pred = c("mean","median"),
                                  var.names = c("y","trial"),
                                  formula = NULL){
  
  point.pred <- match.arg(point.pred)
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  ## How many trials?
  trials <- unique(x[,var.names[2]])
  n.trials <- length(trials)
  ## Need to redefine the individual level
  alpha <- 1 - level
  i.level <- 1 - alpha/neval.1
  ## set up dataset
  x2 <- x[,var.names]
  if(dim(x2)[2] != 2) stop("var.names might be wrong")
  ## do this 'neval' number of times
  pdi.c <- matrix(NA, ncol = 3, nrow = neval.1)
  for(i in 1:neval.1){
    sub.sample <- matrix(nrow = n.trials, ncol = s.size)
    k <- 1
    for(j in trials){
      xs <- x2[x2[,2] == j,var.names[1]]
      sub.sample[k,] <- sample(xs, size = s.size, replace = replace)
      k <- k + 1
    }
    pdi.c[i,] <- pred_int_conformal(as.vector(sub.sample), point.pred = point.pred,
                                    level = level, neval = neval.2)
  }
  ## How should I summarize these sets?
  ## Is this the intersection?
  ## ans <- c(mean(pdi.c[,1]), max(pdi.c[,2]), min(pdi.c[,3]))
  ## Or quantiles?
  ## ans <- c(mean(pdi.c[,1]), quantile(pdi.c[,2], probs = 0.5), quantile(pdi.c[,3],probs = 0.5))
  ## column means?
  ans <- colMeans(pdi.c)
  names(ans) <- c("fit","lwr","upr")
  return(ans)
}