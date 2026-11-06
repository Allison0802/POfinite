\name{pofinite_fit}
\alias{pofinite_fit}
\title{Fit a Fixed Three-Start Finite-Iteration Predictor}
\description{
Fit a complementary log-log mean to raw risk pseudo-observations using a fixed
exact-Hessian trust-region procedure and training-only start selection.
}
\usage{
pofinite_fit(x, y, intercepts)
}
\arguments{
  \item{x}{A finite real numeric training matrix with full column rank and at
  least as many rows as columns. Column names, if supplied, must be unique,
  nonempty and nonmissing.}
  \item{y}{A finite real numeric vector with one raw risk pseudo-observation per
  training row. Values outside \eqn{[0,1]} are allowed and are not clipped.}
  \item{intercepts}{Consecutive leading column indices \code{1:k}. These columns
  must be nonempty, disjoint zero-one cell indicators summing to one in every
  row. At least one additional predictor column is required.}
}
\details{
The objective is \eqn{L(\beta)=\frac{1}{2n}\sum_i(y_i-\mu_i)^2}, with
\eqn{\mu_i=1-\exp\{-\exp(x_i\beta)\}}. This corresponds to constant unit Gaussian
working variance, working independence and unit row weights. It is not Cox
partial likelihood. No inference or cluster-robust covariance is implemented.

Three deterministic starts use the global response mean, the response mean in
each intercept cell, and the cell-mean start plus a prescribed perturbation.
Only start means are clipped to \eqn{[10^{-6},1-10^{-6}]} before taking the link.
Column-RMS scaling and QR give orthogonal training coordinates. The perturbation
has zero intercept coordinates and alternating signs on the remaining
coordinates, with training predictor RMS 0.25 and zero mean in every cell.
The QR pivot must preserve the leading intercept columns, and the original
training matrix must reconstruct its orthogonal coordinates within absolute
error \eqn{10^{-8}}. Numerically ill-conditioned start geometry is rejected
explicitly even if the initial QR rank screen passes.

Each start uses initial trust radius 0.5, maximum radius 8, at most 600 proposals
and 400 accepted steps. Actual and predicted reduction determine acceptance and
radius updates. The original stopping tolerances are fixed; there are no user
optimization controls. Strict-root diagnostics require full derivative rank,
scaled score at most \eqn{10^{-10}}, full-update magnitude at most \eqn{10^{-8}},
full-update prediction change at most \eqn{10^{-9}}, and a Hessian eigenvalue
screen with relative tolerance \eqn{10^{-8}}. They do not gate finite output.

Eligibility requires a controlled termination (\code{STRICT_ROOT},
\code{TRUST_RADIUS_STALL}, \code{TRUST_PROPOSAL_LIMIT}, or
\code{TRUST_ACCEPTED_LIMIT}), no numerical error, finite coefficients,
original-matrix training predictors and risk, finite recomputed loss, and
training coordinate relative error at most \eqn{10^{-8}}. Minimum eligible
training loss selects the endpoint. Losses within
\eqn{10^{-12}\max(1,|L_{min}|)} use priority cell mean, global mean, perturbed
cell mean. Strict convergence and derivative rank are diagnostics only.
New data and held-out outcomes do not participate in selection.

A selected finite endpoint need not solve the estimating equation or attain a
local or global minimum. All endpoints, including failures, remain available.
Rank-deficient, intercept-only, weighted and correlated-working designs are
unsupported. The caller owns design construction, training-only preprocessing,
and alias handling; this function does not build pseudo-observations.
}
\value{
An object of class \code{pofinite} containing:
\item{coefficients}{Named selected coefficients, or missing values if no endpoint is eligible.}
\item{status}{\code{FINITE_ENDPOINT_SELECTED} or \code{NO_ELIGIBLE_ENDPOINT}.}
\item{selected}{Selected start name, or \code{NA_character_}.}
\item{selection}{Selected name, a table of start losses, coordinate errors,
termination reasons, eligibility and selected flags, and eligible loss spread.}
\item{endpoints}{All three solver results, including coefficients, diagnostic
metrics, termination reason, error text, accepted count, proposal trace,
proposals, and accepted coefficient path. A numerical failure may have unavailable diagnostics.}
\item{starts}{The three deterministic starting vectors.}
\item{intercepts, column_names, n_columns, n}{Training design metadata.}
}
\seealso{\code{\link{predict.pofinite}}}
\examples{
z <- seq(-1, 1, length.out = 30)
x <- cbind(cell = 1, z = z)
y <- -expm1(-exp(-0.7 + 0.5 * z))
fit <- pofinite_fit(x, y, intercepts = 1L)
fit
head(predict(fit, x, type = "risk"))
coef(fit)
}
