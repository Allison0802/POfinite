\name{predict.pofinite}
\alias{predict.pofinite}
\alias{coef.pofinite}
\alias{print.pofinite}
\title{Inspect a Finite-Iteration Fit or Predict from Its Selected Endpoint}
\description{Predict using the endpoint selected exclusively from training data,
extract its coefficients, or display the training endpoint dispositions.}
\usage{
\method{predict}{pofinite}(object, newdata,
  type = c("risk", "survival", "link"), ...)
\method{coef}{pofinite}(object, ...)
\method{print}{pofinite}(x, ...)
}
\arguments{
\item{object, x}{A \code{pofinite} fit.}
\item{newdata}{A real numeric matrix with the same columns as the training
matrix. If the training matrix had names, new data must have the same unique
names; their order may differ and is matched automatically. For unnamed
training matrices, new data must also be unnamed and use training column order.}
\item{type}{\code{"risk"} returns \eqn{1-\exp\{-\exp(\eta)\}},
\code{"survival"} returns \eqn{\exp\{-\exp(\eta)\}}, and \code{"link"}
returns \eqn{\eta=x\beta}.}
\item{...}{Reserved for S3 generic compatibility; currently unused.}
}
\details{
Missing (\code{NA} or \code{NaN}) predictor rows return missing predictions;
infinite predictor entries are errors. If finite entries overflow in the
matrix product, the affected predictions are missing. Finite linear predictors
can yield exact zero or one probabilities through floating-point saturation.
No alternate start is selected after prediction failure. Prediction errors do
not change the fitted object. Zero-row new data are supported.

The caller must preserve the fitted design meaning and preprocessing, including
cell indicators and any aliases removed before fitting. Only dimensions and
column naming are checked for new data; extrapolation is not screened. No
standard errors, confidence intervals or robust covariance are available.
}
\value{
\code{predict} returns a numeric vector in input row order, named from new-data
row names when present. It errors if no training endpoint was eligible.
\code{coef} returns the selected coefficient vector (missing if fitting failed).
\code{print} displays status and endpoint selection and returns its input invisibly.
}
\examples{
x <- cbind(cell = 1, z = seq(-1, 1, length.out = 20))
fit <- pofinite_fit(x, rep(0.3, nrow(x)), 1L)
newx <- x[1:3, c("z", "cell")]
newx[2, "z"] <- NA_real_
predict(fit, newx)
predict(fit, newx, type = "survival")
}
