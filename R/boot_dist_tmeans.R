#' Bootstrapped distribution of trial means
#' 
#' @title Trial means bootstrap
#' @param formula standard R formula of the form response ~ trial
#' @param data data frame
#' @param nsim number of bootstrap samples
#' @param var.names alternative to formula
#' @return vector with boostrapped trial means
#' @export
#' @examples 
#' \donotrun{
#' data(soyrs)
#' ## Simply calculate the trial means
#' tmns <- aggregate(lrr ~ Trial_ID, data = soyrs, FUN = mean)
#' ## Bootstrapped stratified trial means
#' btm <- boot_tmeans(lrr ~ Trial_ID, data = soyrs, R = 2e3)$dat
#' 
#' pdi.cf <- pred_int_conformal_df(formula = lrr ~ Trial_ID, x = soyrs)
#' 
#' btm.q <- quantile(btm$ys, probs = c(0.5, 0.025, 0.975))
#' 
#' ggplot() + xlab("lrr") + 
#' geom_density(data = btm, aes(x = ys)) + 
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

boot_tmeans <- function(formula=NULL, data, R = 500, 
                             var.names = c("y","trial"),
                             ci.type = "basic"){
  
  ## Sample with replacement hierarchically
  ## First sample a given trial
  ## Then sample observations within a trial conserving 
  ## Uneven sample sizes
  
  ## First: extract variables
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  
  dat <- data[,var.names]
  names(dat) <- c("y","trial")
  
  tmfb <- function(data, indices){
    ans <- aggregate(formula = y ~ trial, 
                   data = data, 
                   subset = indices, FUN = mean)$y
    ans
  }

  btm <- boot(dat, statistic = tmfb, strata = dat$trial, R = R)
  btm.ci <- boot.ci(btm, type = ci.type)
  
  ans <- list(dat = data.frame(ys = c(btm$t)), boot.ci = btm.ci)
  ans
}




