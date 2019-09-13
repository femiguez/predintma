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

