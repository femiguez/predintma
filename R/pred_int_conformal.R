
#' Calculate prediction interval based on conformal inference. 
#' 
#' @title Prediction interval for a 'numeric' based on conformal inference
#' @name pred_int_conformal
#' @description Conformal inference prediction interval for a numeric
#' @param x should be an object of class "numeric"
#' @param neval number of evaluations to perform
#' @param level coverage level with default 0.95
#' @param point.pred point prediction either 'mean' or 'median'
#' @param method totally experimental either 'quantile' or 'deviation'
#' @return a prediction interval
#' @details I wrote this function after reading this tutorial: https://cdsamii.github.io/cds-demos/conformal/conformal-tutorial.html
#' @export
#' @examples 
#' \dontrun{
#' ## Prediciton interval for simulated data
#' set.seed(123)
#' x <- rnorm(100)
#' pdi.t <- pred_int_tdist(x) 
#' pdi.q <- pred_int_conformal(x)
#' pdi.d <- pred_int_conformal(x, method = "deviation")
#' pdi.j <- pred_int_conformal(x, method = "jackknife")
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
#' Cbounds <- c(u.candidate[match(5, 
#'                                  apply(as.matrix(u.candidate),1,
#'                                  function(x){floor(2.5+100*(mean(c(U, x) < x)))}))],
#'                                  rev(u.candidate)[match(95, rev(apply(as.matrix(u.candidate),1,
#'                                  function(x){ceiling(-2.5+100*(mean(c(U, x) < x)))})))])
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
                               method = c("quantile","deviation",
                                          "jackknife")){
  
  if(class(x) != "numeric" & class(x) != "integer") stop("class should be numeric or integer")
  point.pred <- match.arg(point.pred)
  method <- match.arg(method)
  
  ## Default values in case the methods below fail
  lower.bound <- min(x)
  upper.bound <- max(x)
  
  ## This is our linear evaluation
  xx <- seq(from = min(x), to = max(x), length.out = neval)
  
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
    ## Calculate sum of x
    ## This follows this publication https://arxiv.org/abs/1809.07441
    ## Distribution-free prediciton intervals for random-effects
    ## With an algorithm for small sample sizes
    alpha <- (1 - level)
    pvalues <- numeric(neval)
    length.x <- length(x)
    for(i in 1:neval){
      x.a <- c(x, xx[i]) ## Augmented x
      m.x.a <- mean(x.a) ## mean of augmented x
      ar.xx.i <- abs(xx[i] - m.x.a) ## Residual for each proposed value  
      ar.x <- abs(x - m.x.a) ## Absolute residual for x
      ## Calculate p-value for each u
      pvalues[i] <-  sum(ar.x >= ar.xx.i) / (length.x + 1)
    }
    ## Rule for determining lower bound
    if(pvalues[1] >= alpha/2){
      lower.bound <- xx[1]
    }else{
      ## First half of evaluations
      xx2 <- xx[1:floor(neval/2)]
      pvalues2 <- pvalues[1:floor(neval/2)]
      lower.bound <- max(xx2[pvalues2 < alpha/2])
    }
    ## Rule for determining upper bound
    if(pvalues[neval] >= alpha/2){
      upper.bound <- xx[neval]
    }else{
      xx3 <- xx[floor(neval/2):neval]
      pvalues3 <- pvalues[floor(neval/2):neval]
      upper.bound <- min(xx3[pvalues3 < alpha/2])
    }
  }
  
  if(method == "jackknife"){
    rs.x <- numeric(length(x))
    for(i in 1:length(x)){
      x.m1 <- x[-i]
      m.x.m1 <- mean(x.m1)
      rs.x[i] <- abs(x[i] - m.x.m1)
    }
    k <- ceiling(length(x)*level)
    dd <- sort(rs.x)[k]
    ## print(rs.x)
    ## print(dd)
    lower.bound <- mean(x) - dd
    upper.bound <- mean(x) + dd
  }
  
  if(point.pred == "median"){
    pprd <- median(x)
  }else{
    pprd <- mean(x)
  }
  ans <- c(pprd, lower.bound, upper.bound)
  names(ans) <- c("fit","lwr","upr")
  return(ans)
}