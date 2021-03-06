% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_profile.R
\name{pred_profile}
\alias{pred_profile}
\title{pred_profile}
\usage{
pred_profile(
  x,
  probs = seq(0.01, 0.99, by = 0.01),
  method = c("tdist", "conformal", "quantile"),
  m.method = c("quantile", "deviation", "jackknife"),
  neval = 200,
  point.pred = c("mean", "median")
)
}
\arguments{
\item{probs}{probabilities at which to calculate a probability interval}

\item{method}{either "tdist", "conformal" or "quantile"}

\item{m.method}{conformal method (quantile, deviation or jackknife)}
}
\value{
an object of class 'pred_profile'
}
\description{
Prediction Interval Profiles
}
\examples{
\dontrun{
set.seed(12345)
x <- rnorm(25)
## t-dist method
pp1 <- pred_profile(x)
## conformal - quantile
pp2 <- pred_profile(x, method = "conformal")
## conformal - deviation
pp3 <- pred_profile(x, method = "conformal", m.method = "deviation")
## conformal - jackknife
pp4 <- pred_profile(x, method = "conformal", m.method = "jackknife")
par(mfrow = c(2,2))
plot(pp1)
plot(pp2)
plot(pp3)
plot(pp4)
par(mfrow=c(1,1))
}

}
