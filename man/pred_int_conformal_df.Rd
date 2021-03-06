% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_int_conformal_df.R
\name{pred_int_conformal_df}
\alias{pred_int_conformal_df}
\title{Conformal Prediction Interval for Random-effects Meta-analysis}
\usage{
pred_int_conformal_df(
  x,
  method = "subsampling",
  level = 0.95,
  m.method = c("quantile", "deviation", "jackknife"),
  neval.1 = 500,
  neval.2 = 500,
  s.size = 1,
  replace = FALSE,
  point.pred = c("mean", "median"),
  var.names = c("y", "trial"),
  formula = NULL
)
}
\arguments{
\item{x}{should be an object of class "data.frame"}

\item{method}{"subsampling"}

\item{level}{coverage level with default 0.95}

\item{m.method}{conformal method used either: 'quantile','deviation' or 'jackknife'}

\item{neval.1}{number of evaluations at the 'trial' level.}

\item{neval.2}{number of evaluations at the 'rep' level.}

\item{s.size}{within-group sample size, when the number of groups is less than 20, it should be higher than 1.}

\item{replace}{replacement option in the within-group sampling}

\item{point.pred}{'point prediction', either 'mean' or 'median'}

\item{var.names}{variable names to be passed to the 'data.frame' methods}

\item{formula}{alternative formula interface, if this is used, 'var.names' is ignored}
}
\value{
a prediction interval for a "new_trial"
}
\description{
Conformal Prediction for Random-effects Model using Subsampling
}
\details{
The idea for this method comes from the paper "Distribution-free Prediction Sets with Random Effects" by Dunn and Wasserman https://arxiv.org/abs/1809.07441
}
\examples{
\dontrun{
## Using soybean row spacing
data(soyrs)
pdi.conf <- pred_int_conformal_df(formula = lrr ~ Trial_ID, x = soyrs) 
}



}
