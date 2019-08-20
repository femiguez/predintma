#' Soybean row spacing on-farm trials
#'
#' A dataset containing the Year, trial id, yield of treatment, yield of control,
#' logarithm of the response ratio and the 'id'
#'
#' @format A data frame with 119 rows and 6 variables:
#' \describe{
#'   \item{Year}{ -integer- year that the trial was conducted.}
#'   \item{Trial_ID}{ -factor- given name for the specific trial.}
#'   \item{TRT1_Yld}{ -numeric- yield for 'narrow' (15 in) row spacing.}
#'   \item{TRT2_Yld}{ -numeric- yield for the 'wide' (30 in) row spacing.}
#'   \item{lrr}{ -numeric- logarithm of the 'response ratio' i.e. log(TRT1_Yld/TRT2_Yld).}
#'   \item{id}{ -integer- added 'id' column which orders observations by the magnitude of
#'              the 'lrr'}
#' }
#' @source \url{https://analytics.iasoybeans.com/cool-apps/ISOFAST/}
"soyrs"

 #' Import packages needed for predintma to work correctly
 #' @import coda emmeans knitr lme4 metafor MCMCglmm stats
 NULL