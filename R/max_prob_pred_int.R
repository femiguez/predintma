#' Calculate maximum probability for the prediction of a small sample size using conformal prediction
#' 
#' 
#' @title maximum probability for the prediciton interval 'small sample size'
#' @name max_prob_pred_int2
#' @description Calculates the maximum probability for a prediction interval 
#' @param x should be a vector
#' @param n number of evaluations to perform in the line search of probabilities
#' @param neval number of evaluations in the conformal algorithm (see pred_int_conformal)
#' @param tol tolerance default is (0.001) see details
#' @param method either 'tdist' (assumes normality) or 'conformal' (distribution-free)
#' @param m.method method used to compute conformal prediction interval: either "quantile", "deviation" or "jackknife"
#' @return a single value which represents a 'suggestion' for the maximum level of probability given the data
#' @details The idea is to find the maximum level of probability that will produce a conformal prediction interval
#'          which matches the minimum and maximum values in the observed sample. This is approximate and it depends on 
#'          the tolerance. The distance is calculated as abs(x.min - calc.lower.bound) + abs(x.max - calc.upper.bound).
#'          The deviation tolerance is calculated as a proportion of the range in the observed sample. i.e. (x.max - x.min)*tol
#' @note At this moment, with current testing, for method "tdist" this is not reliable and "max_prob_pred_int" 
#'       should be used. 
#' @export
#' @examples 
#' \dontrun{
#' set.seed(12345)
#' mp.pdi <- max_prob_pred_int2(rnorm(10), method = "conformal") 
#' ## For this particular sample, the 'suggested' maximum is about 82%
#' }
#'
#'

max_prob_pred_int2 <- function(x, n = 200, neval = 200, tol = 0.001, 
                              method = c("tdist","conformal"),
                              m.method = c("quantile", "deviation","jackknife")){
  
  ## Not suitable for a very small sample size
  if(length(x) < 5) stop("sample size should be at least 5")
  
  method <- match.arg(method)
  m.method <- match.arg(m.method)
  
  prob.level <- seq(0.001,1,length.out = n)
  x.max <- max(x)
  x.min <- min(x)
  
  ans <- 0
  
  for(i in 1:n){
    if(method == "tdist"){
      pdi <- pred_int(x, level = prob.level[i])
    }else{
      pdi <- pred_int_conformal(x, neval = neval, method = m.method, level = prob.level[i])
    }
    dv <- abs(x.min - pdi[2]) + abs(x.max - pdi[3])
    if(dv < tol){
      ans <- prob.level[i]
      break
    }
  }
  return(ans)
}


#' Optimize the maximum probability for the prediction of a small sample size using either the 
#' t-distribution or conformal methods
#' 
#' @title Maximum probability for the prediciton interval 'small sample size'
#' @name max_prob_pred_int
#' @description Finds the maximum probability for a prediction interval using optimization 
#' @param x should be a vector
#' @param method either 'tdist' (assumes normality), 'conformal' (distribution-free), or non-parametric ('npar')
#' @param m.method method used to compute conformal prediction interval: either "quantile", "deviation" or "jackknife"
#' @param interval maximum and minimum values for the optimization search
#' @param alpha.penalty whether to include a penalty for alpha
#' @param scale whether to scale the input vector
#' @return a single value which represents a 'suggestion' for the maximum level of probability given the data
#' @details The idea is to find the maximum level of probability that will produce a prediction interval
#'          which matches the minimum and maximum values in the observed sample.  The distance is calculated as abs(x.min - calc.lower.bound) + abs(x.max - calc.upper.bound).
#'          
#' @export
#' @examples 
#' \dontrun{
#' set.seed(12345)
#' x <- rnorm(10)
#' mp.pdi <- max_prob_pred_int(x) 
#' ## For this particular sample, the maximum is about 76%
#' data(soyrs)
#' soyrs.s <- aggregate(lrr ~ Trial_ID, data = soyrs, FUN = mean)
#' mp.pdi.soy.tdist <- max_prob_pred_int(soyrs.s$lrr)
#' mp.pdi.soy.conf <- max_prob_pred_int(soyrs.s$lrr, method = "conformal")
#' print(round(c(mp.pdi.soy.tdist, mp.pdi.soy.conf1),2))
#' }
#'
#'
max_prob_pred_int <- function(x, method = c("tdist","conformal","npar"),
                                  m.method = c("quantile","deviation","jackknife"),
                                  interval = c(0,1), alpha.penalty = 0,
                                  scale = FALSE){
  
  method <- match.arg(method)
  m.method <- match.arg(m.method)
  
  if(method == "tdist"){
    ans <- optimize(f = mpdi_obj, interval = interval, 
                    x = x, method = method, m.method = m.method)$minimum
  }
  
  if(method == "conformal"){
    ## It turns out that for the conformal method this optimization
    ## problem is harder, but not that hard
    ## First line-search the optimization function
    alphas <- seq(0.05, 0.95, by = 0.05)
    objf <- numeric(length(alphas))
    for(i in 1:length(objf)){
      objf[i] <- mpdi_obj(alphas[i], x = x, 
                          method = method, m.method = m.method)
    }
    if(abs(max(objf) - min(objf)) < 0.001){
      stop("objective function is flat. Can't optimize it.")
    } 
    alob <- data.frame(alpha = alphas, objfv = objf)
    mm.alpha <- max(alob[alob$objfv == min(objf),"alpha"])
    ans <- optimize(f = mpdi_obj, interval = c(mm.alpha-0.04,mm.alpha+0.04), 
                    x = x, method = method, m.method = m.method,
                    alpha.penalty = alpha.penalty,
                    scale = scale)$minimum
  }
  
  if(method == "npar"){
    ans <- 1 - (length(x) - 1)/(length(x) + 1)
  }
  
  return(1 - ans)
}

#' Objective funciton for optimization
#' 
#' @title objective function for probability of a prediciton interval for 'small' sample sizes.
#' @name mpdi_obj
#' @description objective function for probability of a prediciton interval for 'small' sample sizes. 
#' @param alpha miscoverage or 'error rate'
#' @param x a vector
#' @param method either 'tdist' (assumes normality) or 'conformal' (distribution-free)
#' @param m.method method used to compute conformal prediction interval: either "quantile", "deviation" or "jackknife"
#' @param alpha.penalty whether to include an alpha penalty (default 0 or 'no')
#' @param scale whether to scale the input vector. This only makes sense if the alpha.penalty is different from zero. 
#' @return a single value which represents a value that should be minimized
#' @details The idea is to find the maximum level of probability that will produce a prediction interval
#'          which matches the minimum and maximum values in the observed sample.  The distance is calculated as abs(x.min - calc.lower.bound) + abs(x.max - calc.upper.bound).
#'          
#' @export
#' @examples 
#' \dontrun{
#' set.seed(12345)
#' x <- rnorm(10)
#' alphas <- seq(0,1, 0.05)
#' objf <- numeric(length(alphas))
#' for(i in 1:length(objf)){
#'   objf[i] <- mpdi_obj(alphas[i], x = x, method = "conformal")
#' }
#' qplot(alphas, objf, geom = "line")
#' ## Trying the t-distribution
#' y <- rt(10, df = 1)
#' alphas <- seq(0,1, 0.05)
#' objf <- numeric(length(alphas))
#' for(i in 1:length(objf)){
#'   objf[i] <- mpdi_obj(alphas[i], x = y, method = "conformal")
#' }
#' qplot(alphas, objf, geom = "line")
#' }

mpdi_obj <- function(alpha, x, method = c("tdist","conformal"),
                     m.method = c("quantile","deviation","jackknife"),
                     alpha.penalty=0, scale = FALSE){
  
  method <- match.arg(method)
  m.method <- match.arg(m.method)
  ## This function will return the objective 
  ## for the minimization
  ## Component 1 is simply alpha
  level <- 1 - alpha
  if(scale) x <- scale(x)[,1]
  ## Component two is the deviation
  x.max <- max(x)
  x.min <- min(x)
  if(method == "tdist"){
    pdi <- pred_int(x, level = level)
    ans <- abs(x.max - pdi[3]) + abs(x.min - pdi[2])
  }
  
  if(method == "conformal"){
    pdi <- pred_int_conformal(x, level = level, method = m.method)
    ans <- abs(x.max - pdi[3]) + abs(x.min - pdi[2]) - alpha * alpha.penalty 
  }
  ## name the result
  setNames(ans, "obj_fun")
  return(ans)
}
