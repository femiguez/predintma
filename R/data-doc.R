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

#' Soybean Stratego Foliar Fungicide
#'
#' A dataset containing the Year, trial id, rep, yield of treatment, yield of control and the
#' logarithm of the response ratio 
#'
#' @format A data frame with 200 rows and 6 variables:
#' \describe{
#'   \item{Year}{ -integer- year that the trial was conducted.}
#'   \item{Trial_ID}{ -factor- given name for the specific trial.}
#'   \item{Rep}{ -integer- rep number.}
#'   \item{TRT_Yld}{ -numeric- yield for treated.}
#'   \item{CTR_Yld}{ -numeric- yield for untreated (control).}
#'   \item{lrr}{ -numeric- logarithm of the 'response ratio' i.e. log(TRT_Yld/CTR_Yld).}
#' }
#' @source \url{https://analytics.iasoybeans.com/cool-apps/ISOFAST/}
"soysff"

#' Maize Foliar Fungicide
#'
#' A dataset containing the trial id, year, rep, yield of treatment, yield of control and the
#' logarithm of the response ratio 
#'
#' @format A data frame with 703 rows and 6 variables:
#' \describe{
#'   \item{Trial_ID}{ -factor- given name for the specific trial.}
#'   \item{Year}{ -integer- year that the trial was conducted.}
#'   \item{Rep}{ -integer- rep number.}
#'   \item{TRT_Yld}{ -numeric- yield for treated (fungicide Headline).}
#'   \item{CTR_Yld}{ -numeric- yield for untreated (control).}
#'   \item{lrr}{ -numeric- logarithm of the 'response ratio' i.e. log(TRT_Yld/CTR_Yld).}
#' }
#' @source \url{https://analytics.iasoybeans.com/cool-apps/ISOFAST/}
"mzfg"

#' Maize Foliar Fungicide (Stratego)
#'
#' A dataset containing the trial id, year, rep, yield of treatment, yield of control and the
#' logarithm of the response ratio 
#'
#' @format A data frame with 153 rows and 4 variables:
#' \describe{
#'   \item{Trial_ID}{ -factor- given name for the specific trial.}
#'   \item{TRT_Yld}{ -numeric- yield for treated (fungicide Stratego).}
#'   \item{CTR_Yld}{ -numeric- yield for untreated (control).}
#'   \item{lrr}{ -numeric- logarithm of the 'response ratio' i.e. log(TRT_Yld/CTR_Yld).}
#' }
#' @source \url{https://analytics.iasoybeans.com/cool-apps/ISOFAST/}
"mzstr"

#' Soybean Fungicide (Priaxor)
#'
#' A dataset containing the trial id, year, rep, yield of treatment, yield of control and the
#' logarithm of the response ratio 
#'
#' @format A data frame with 191 rows and 6 variables:
#' \describe{
#'   \item{Trial_ID}{ -factor- given name for the specific trial.}
#'   \item{Year}{ -integer- year when the trial was conducted.}
#'   \item{Rep}{-integer- replication.}
#'   \item{TRT_Yld}{ -numeric- yield for treated (fungicide Priaxor).}
#'   \item{CTR_Yld}{ -numeric- yield for untreated (control).}
#'   \item{lrr}{ -numeric- logarithm of the 'response ratio' i.e. log(TRT_Yld/CTR_Yld).}
#' }
#' @source \url{https://analytics.iasoybeans.com/cool-apps/ISOFAST/}
"soyprx"

 #' Import packages needed for predintma to work correctly
 #' @import boot coda emmeans knitr lme4 metafor MCMCglmm parallel stats
 NULL