#' Prediction Interval Profiles
#' 
#' @title pred_profile
#' @param probs probabilities at which to calculate a probability interval
#' @param method either "tdist" or "conformal"
#' @param m.method conformal method (quantile, deviation, jackknife)
#' @return data frame with prediction intervals at specified probs
#' @export
#' @examples 
#' \donotrun{
#' set.seed(12345)
#' x <- rnorm(25)
#' ## t-dist method
#' pp1 <- pred_profile(x)
#' ## conformal - quantile
#' pp2 <- pred_profile(x, method = "conformal")
#' ## conformal - deviation
#' pp3 <- pred_profile(x, method = "conformal", m.method = "deviation")
#' ## conformal - jackknife
#' pp4 <- pred_profile(x, method = "conformal", m.method = "jackknife")
#' par(mfrow = c(2,2))
#' plot(pp1)
#' plot(pp2)
#' plot(pp3)
#' plot(pp4)
#' par(mfrow=c(1,1))
#' }
#' 

pred_profile <- function(x, probs = seq(0.01, 0.99, by = 0.01),
                         method = c("tdist","conformal","quantile"),
                         m.method = c("quantile","deviation","jackknife"),
                         neval = 200, point.pred = c("mean","median")){
  
  method <- match.arg(method)
  m.method <- match.arg(m.method)
  
  ans <- data.frame(probs = probs, m = NA, lb = NA, ub = NA)
  
  l.probs <- length(probs)
  
  for(i in 1:l.probs){
    
    if(method == "tdist"){
      ans[i,c("m","lb","ub")] <- pred_int(x, level = probs[i])
    }
    
    if(method == "conformal"){
        ans[i,c("m","lb","ub")] <- pred_int_conformal(x, neval = neval, level = probs[i],
                                                      point.pred = point.pred, method = m.method)
    }
    
    if(method == "quantile"){
      alpha <- 1 - probs[i]
      half.alpha <- alpha/2
      ans[i,c("m","lb","ub")] <- quantile(x, probs = c(0.5, half.alpha, c(1 - half.alpha)))
    }
    
  }
  ret <- list(ans = ans, method = method, m.method = m.method, data = x)
  class(ret) <- c("pred_profile")
  ret
}

#' Plot Prediction Interval Profiles
#' 
#' @title plot.pred_profile
#' @param max.prob logical whether to plot the optmiized maximum probability
#' @param main optional title for the plot
#' @param max.prob.col color for the plotting of the maximum probability
#' @return profile plots
#' @export
#' 
plot.pred_profile <- function(x, max.prob = TRUE, main = NULL, max.prob.col = "blue"){
  
  if(missing(main)){
    if(x$method != "conformal") ptitle <- x$method
    if(x$method == "conformal") ptitle <- paste(x$method, x$m.method)
  }

  if(max.prob == TRUE & x$method == "quantile"){
    mpp <- opt_max_prob_pred_int(x$data, method = "npar")
  } 
  
  if(max.prob == TRUE & x$method != "quantile"){
    mpp <- opt_max_prob_pred_int(x$data, method = x$method, m.method = x$m.method)
  }
  
  tmp <- x$ans
  
  plot(x = tmp$m, y = tmp$probs, type = "p", pch = ".",
       xlim = c(min(tmp$lb),max(tmp$ub)),
       main = ptitle,
       xlab = "y", ylab = "probabilities")
  lines(x = tmp$lb, y = tmp$probs)
  lines(x = tmp$ub, y = tmp$probs)
  abline(h = mpp, col = max.prob.col)
  
}



