#' Bootstrapped distribution of trial means
#' 
#' @title Bootstrap distribution of trial means
#' @name boot_tmeans
#' @param formula standard R formula of the form response ~ trial
#' @param data data frame
#' @param R number of bootstrap samples
#' @param var.names alternative to formula
#' @param ci.type type of confidence interval (see boot::boot.ci)
#' @param level prediction level
#' @param ncpus number of cpus to use, passed to boot, negligible benefit for small datasets
#' @return list with 'dat': data frame with the bootstrapped means, \cr 
#'          'boot.ci': confidence interval computed by 'boot.ci' \cr
#'          and 'pdi': prediction interval based on quantiles 
#' @export
#' @examples 
#' \dontrun{
#' require(ggplot2)
#' data(soyrs)
#' ## Simply calculate the trial means
#' tmns <- aggregate(lrr ~ Trial_ID, data = soyrs, FUN = mean)
#' ## Bootstrapped stratified trial means
#' btm <- boot_tmeans(lrr ~ Trial_ID, data = soyrs, R = 2e3)
#' 
#' pdi.cf <- pred_int_conformal_df(formula = lrr ~ Trial_ID, x = soyrs)
#' 
#' btmd <- btm$dat
#' btm.q <- btm$pdi
#' 
#' ggplot() + xlab("lrr") + 
#' geom_density(data = btmd, aes(x = ys)) + 
#' geom_jitter(data = soyrs, aes(x = lrr, y = 2.5)) + 
#' geom_jitter(data = tmns, aes(x = lrr, y = 5), color = "blue", size = 1.2) + 
#' geom_point(aes(x = pdi.cf[1], y = -1), color = "orange", size = 1.2) + 
#' geom_errorbarh(mapping = aes(xmin = pdi.cf[2], xmax = pdi.cf[3],
#'                               y = -1), color = "orange", size = 1.2) +
#' geom_point(aes(x = btm.q[1], y = -2), color = "red", size = 1.2) + 
#' geom_errorbarh(aes(xmin = btm.q[2], xmax = btm.q[3], y = -2),
#'                color = "red", size = 1.2)
#' 
#' }

boot_tmeans <- function(formula = NULL, data, R = 500, 
                        var.names = c("y","trial"),
                        ci.type = "basic", 
                        level = 0.95,
                        ncpus = 1){
  
  ## Sample with replacement with strata
  alpha <- 1 - level
  half.alpha <- alpha/2
  
  if(ncpus > 1){
    prll <- "multicore"
  }else{
    prll <- "no"
  }
  
  ## First: extract variables
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))

  dat <- data[,var.names]
  names(dat) <- c("y","trial")
  dat$trial <- as.factor(dat$trial)

  tmfb <- function(data, indices){
    ans <- aggregate(formula = y ~ trial, 
                     data = data,
                     subset = indices,
                     FUN = mean)$y
    ans
  }

  btm <- boot::boot(dat, statistic = tmfb, strata = dat$trial, R = R,
                    sim = "ordinary", stype = "i",
                    parallel = prll, ncpus = ncpus)
  
  btm.ci <- boot::boot.ci(btm, conf = level, type = ci.type)
  
  pdi <- quantile(c(btm$t), probs = c(0.5, half.alpha, 1 - half.alpha))
  
  ans <- list(dat = data.frame(ys = c(btm$t)), 
              boot.ci = btm.ci,
              pdi = pdi)
  ans
}


#' Bootstrap-based prediction intervals 
#' 
#' @title Bootstrap-based prediction intervals
#' @name pred_int_boot
#' @param formula standard R formula of the form response ~ trial
#' @param data data frame
#' @param level prediction level
#' @param R number of bootstrap samples
#' @param ncpus number of cpus to use, passed to boot, negligible benefit for small datasets
#' @return prediction interval
#' @details a simple wrapper for boot_tmeans
#' @export
#' 

pred_int_boot <- function(formula, data, level = 0.95, R = 500, ncpus=1, ...){
  pdi <- boot_tmeans(formula = formula, data = data, R = R, 
                     level = level, ncpus = ncpus, ...)$pdi
  pdi
}

#' This will implement hierarchical bootstrapping as oppossed to
#' stratified bootstrapping, which is the default
#' 
#' @title hierarchical bootstrap
#' @name boot_tmeans2 
#' @param formula standard R formula, as in resp ~ trial
#' @param data data frame
#' @param R number of replicates
#' @param var.names optional alternative to formula
#' @param level probability level
#' @param ncpus number of cpus for parallel resampling (not implemented yet)
#' @export
#' @examples 
#' \dontrun{
#' require(ggplot2)
#' ## How does it compare to stratified resampling?
#' data(soyrs)
#' btm1 <- boot_tmeans(lrr ~ Trial_ID, data = soyrs)
#' btm2 <- boot_tmeans2(lrr ~ Trial_ID, data = soyrs)
#' 
#' dat <- data.frame(method = rep(c("boot","boot2"), each = nrow(btm2$dat)),
#'                   pos = rep(c(2,5), each = nrow(btm2$dat)),
#'                   ys = c(btm1$ys, btm2$dat$ys))
#'                   
#' ggplot(data = dat) + 
#'      xlab("lrr") +
#'      geom_point(aes(x = ys, y = pos, color = method)) + 
#'      geom_jitter(aes(x = ys, y = pos, color = method)) +
#'      geom_density(aes(x = ys, color = method))
#'      
#' ## The two methods give nearly identical answers
#' 
#' }
#' 

boot_tmeans2 <- function(formula=NULL, data, R = 500, 
                         var.names = c("y","trial"),
                         level = 0.95, ncpus = 1){
  
  alpha <- 1 - level
  half.alpha <- alpha/2
  
  if(ncpus > 1){
    prll <- "multicore"
  }else{
    prll <- "no"
  }
  
  ## First: extract variables
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  
  dat <- data[,var.names]
  names(dat) <- c("y","trial")
  
  trials <- unique(dat$trial)
  n.trials <- length(trials)
  
  ## Bootstrapped trial mean
  btm <- function(x, i){
    ans <- mean(x[i])
  }
  
  res <- data.frame(iter = rep(1:R, each = n.trials), ys = NA)
  
  for(i in 1:R){
    s.trials <- sample(trials, size = n.trials, replace = TRUE)
    tmp.m <- numeric(n.trials)
    
    for(j in 1:n.trials){
      tmp <- subset(dat, trial == s.trials[j])$y
      tmp.m[j] <- boot(tmp, btm, R = 1)$t
    }
    
    i1 <- (i - 1)*n.trials + 1
    i2 <- i * n.trials
    res[i1:i2,"ys"] <- tmp.m    
  }
  
  pdi <- quantile(res$ys, probs = c(0.5, half.alpha, 1 - half.alpha))
  
  ans <- list(dat = res, pdi = pdi)
  return(ans)
}


#' Bootstrap-based prediction intervals (method 2)
#' 
#' @title Bootstrap-based prediction intervals (method 2)
#' @name pred_int_boot2
#' @param formula standard R formula of the form response ~ trial
#' @param data data frame
#' @param level prediction level
#' @param R number of bootstrap samples
#' @param ncpus number of cpus to use, passed to boot, negligible benefit for small datasets
#' @return prediction interval
#' @details a simple wrapper for boot_tmeans2
#' @export
#' 
pred_int_boot2 <- function(formula, data, level = 0.95, R = 500, ncpus = 1){
  
  ans <- boot_tmeans2(formula = formula, data = data, 
                      level = level, R = R, ncpus = ncpus)$pdi
  ans
}


