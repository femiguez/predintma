#' Custom function to summarize trial data 
#' 
#' 
#' @title Summarize a data.frame according to a given 'trial'
#' @name tsum
#' @description Summarize data according to a 'trial' variable using either individual
#'  confidence intervals, model-based confidence intervals ('lm') or 'min-max'
#' @param x should be an object of class "data.frame"
#' @param var.names names of variables 'y' (response) and 'trial' variable to summarize by.
#' @param method 'ici' individual confidence interval, 'lm' model-based confidence interval 
#' (with typical assumptions of a linear model), 'min-max' minimum and maximum values
#' @param level coverage level with default 0.95
#' @param order whether to order the summary by mean of response
#' @param add.id whether to add the 'id' label
#' @param formula alternative formula interface, if this is used 'var.names' is ignored.
#' @return a data.frame with the indicated summaries
#' @details this function is for somewhat of narrow use case, but similar things can be accomplished
#' using 'tidyverse' type operations. It uses 'aggregate'.
#' @export
#' @examples 
#' \dontrun{
#' ## Summarize trial data
#' data(soyrs)
#' soyrs.s <- tsum(soyrs, var.names = c("lrr","Trial_ID")) 
#' 
#' ggplot(data = soyrs.s, aes(x = m, y = id)) + 
#' geom_point() + xlab("log RR") + ylab("ID") +
#' geom_vline(xintercept = 0)
#' 
#' }

tsum <- function(x, var.names = c("y","trial"), 
                 method = c("ici","lm","min-max"),
                 level = 0.95, order = TRUE,
                 add.id = TRUE,
                 formula = NULL){
  
  if(class(x) != "data.frame") stop("only for data.frames")
  
  method <- match.arg(method)
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  
  if(method == "ici"){
    ## Calcualte mean by trial
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    ans <- aggregate(formula = frm, data = x, FUN = mean)
    ans$lb <- aggregate(formula = frm, data = x, FUN = lb_fun, level = level)[,2]
    ans$ub <- aggregate(formula = frm, data = x, FUN = ub_fun, level = level)[,2]
    ans$n <- aggregate(formula = frm, data = x, FUN = length)[,2]
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "lm"){
    ## Calcualte intervals using lm
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    nt <- aggregate(formula = frm, data = x, FUN = length)[,2]
    fit <- lm(formula = frm, data = x)
    ndat <- data.frame(unique(x[,var.names[2]]))
    names(ndat) <- var.names[2]
    prd <- predict(fit, newdata = ndat, interval = "confidence")
    ans <- cbind(ndat, prd)
    ans$n <- nt
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "min-max"){
    ## Calcualte mean by trial
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    ans <- aggregate(formula = frm, data = x, FUN = mean)
    ans$lb <- aggregate(formula = frm, data = x, FUN = min)[,2]
    ans$ub <- aggregate(formula = frm, data = x, FUN = max)[,2]
    ans$n <- aggregate(formula = frm, data = x, FUN = length)[,2]
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(order & add.id){
    ans2 <- ans[order(ans[,"m"]),]
    ans2$id <- 1:nrow(ans2)
    return(ans2)
  }else{
    return(ans)
  }
}

lb_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  lbv <- m.x -  tv * se.x
  return(lbv)
}

ub_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  ubv <- m.x +  tv * se.x
  return(ubv)
}

#' @title Proportion of values in 'x' which are included in and interval
#' @name includes
#' @description Count how many elements of vector 'x' are included in a given interval.
#' @param x should be an object of class "numeric"
#' @param lb lower bound
#' @param ub upper bound
#' @return a single value indicating the proportion of elements in 'x' which are contained in
#' a given interval [lb,ub]
#' @details if an element is contained in the interval [lb,ub] it will be counted.
#' https://www.varsitytutors.com/hotmath/hotmath_help/topics/interval-notation
#' @export
#' @examples 
#' \dontrun{
#' x <- rnorm(100)
#' includes(x, -1.8, 1.8)
#' }

includes <- function(x, lb, ub){
  
  if(ub < lb) stop("upper bound should be greater than lower bound")
  tmp <- 0
  for(i in 1:length(x)){
    if(x[i] >= lb & x[i] <= ub){
      tmp <- tmp + 1
    }
  }
  ans <- tmp/length(x)
  ans
}