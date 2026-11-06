# POfinite — finite-iteration pseudo-observation prediction

A fixed three-start finite-iteration predictor for raw risk pseudo-observations.
The matrix API uses a complementary log-log conditional mean, constant unit
Gaussian working variance, working independence, and unit row weights. It has
no non-base runtime dependencies beyond the recommended `stats` and `utils`
packages shipped with R.

Author and maintainer: Yumei Yang.

Version 0.1.0 is an experimental release. It provides a safeguarded
finite-iteration predictor and numerical diagnostics. It does not guarantee
convergence, replace general GEE software, or provide validated coefficient
inference. All supplied examples and tests use synthetic inputs.

Released under the [MIT license](LICENSE.md). Maintainer contact:
yumeiy@email.unc.edu. This package is distributed through GitHub and is not on CRAN.

## Installation

Install the versioned GitHub release:

```r
install.packages("remotes")  # only if needed
remotes::install_github("Allison0802/POfinite@v0.1.0")
```

Alternatively, download `POfinite_0.1.0.tar.gz` from the
[release page](https://github.com/Allison0802/POfinite/releases/tag/v0.1.0)
and install it locally:

```r
install.packages("POfinite_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Why retain finite endpoints?

POfinite is designed to retain finite coefficient estimates and predictions even
when strict-root convergence checks fail. This can be useful when a conventional
GEE fitting workflow reports numerical or convergence difficulties. It evaluates
three fixed starts, selects an eligible finite endpoint by training loss, and
retains the failure and convergence diagnostics for inspection. If no endpoint
is eligible, it reports failure rather than returning substitute predictions.

In a synthetic comparison with geepack 1.3.13 (R 4.4.0), all three all-zero
boundary fits from geepack reported native `error = 1`; POfinite retained finite
coefficients and predictions while also failing its strict-root checks. Both
packages returned finite values in that comparison. Across all five synthetic
fixtures, all 15 endpoints from each package were finite. These results support
POfinite's explicit handling of finite stopped endpoints, not a claim that
geepack necessarily returns nonfinite values or that POfinite always succeeds.

The following synthetic boundary example exercises this behavior without any
clinical data:

```r
library(POfinite)
x0 <- cbind(cell = 1, z = seq(-1.8, 1.8, length.out = 60))
f0 <- pofinite_fit(x0, rep(0, nrow(x0)), intercepts = 1L)
f0$status
f0$selection$table
vapply(f0$endpoints, function(e) e$metrics[["strict_pass"]], numeric(1))
stopifnot(all(is.finite(coef(f0))), all(is.finite(predict(f0, x0))))
```

For an optional geepack comparison using the same starts and unit Gaussian
working variance (geepack must be installed separately):

```r
lapply(f0$starts, function(b0) {
  g <- geepack::geese.fit(
    x0, rep(0, nrow(x0)), id = rep(seq_len(20), each = 3),
    b = b0, gm = 1, family = gaussian(), mean.link = "cloglog",
    variance = "gaussian", corstr = "independence", scale.fix = TRUE,
    scale.value = 1,
    control = geepack::geese.control(epsilon = 1e-8, maxit = 100)
  )
  list(native_error = g$error, coefficients = g$beta,
       finite = all(is.finite(g$beta)))
})
```

Results can vary with software versions and numerical environment. Finite
coefficients or predictions do not establish convergence, a minimum, predictive
accuracy, or valid inference. POfinite is not a general replacement for geepack
and does not provide robust standard errors or confidence intervals.

## Synthetic example

```r
library(POfinite)
z <- seq(-1, 1, length.out = 30)
x <- cbind(cell = 1, z = z)
y <- -expm1(-exp(-0.7 + 0.5 * z))
fit <- pofinite_fit(x, y, intercepts = 1L)
fit
predict(fit, x, type = "risk")
predict(fit, x, type = "survival")
coef(fit)
fit$selection$table
```

`x` must be a finite real numeric matrix. `y` is a finite numeric vector of
unmodified raw risk pseudo-observations; values outside [0,1] are allowed.
The function does not construct pseudo-observations, patient folds or design
matrices. The caller must preserve the intended prediction target and apply
training-only preprocessing.

`intercepts` identifies the first k columns: disjoint 0/1 cell indicators
summing to one in each row, with every training cell represented. The retained
matrix must have full column rank, n >= p, and p > k. The frozen perturbation
is undefined for intercept-only designs, which are explicitly unsupported.
Exact aliases must be removed before fitting. Zero or overflowing column scales
are rejected. QR must preserve the leading intercept block and reconstruct the
original design coordinates within absolute error 1e-8. Numerically unstable
start geometry is an explicit unsupported-design error even when the initial
rank screen passes. No automatic reparameterization changes the prescribed rule.

## Fixed algorithm and endpoint selection

The objective is half the mean squared difference between raw responses and
`-expm1(-exp(x %*% beta))`. Start means alone are clipped to [1e-6, 1-1e-6]
before applying the complementary log-log link. Three deterministic starts are:

1. The global raw-response mean.
2. Each intercept cell's raw-response mean.
3. The cell-mean start plus an orthogonal perturbation with training predictor
   RMS 0.25 and zero cell means. Its non-intercept orthogonal coordinates
   alternate positive and negative signs.

The solver uses the exact objective Hessian in column-RMS-scaled QR coordinates.
Its trust radius begins at 0.5 and is capped at 8. There are at most 600 proposals
and 400 accepted steps per start. A step is accepted only for positive actual
reduction and an actual/predicted reduction ratio above 0.1. A ratio below 0.25
(or nonfinite ratio) divides the radius by four; a ratio above 0.75 with step
norm at least 90% of the radius doubles it up to the cap. The radius stall
threshold is 1e-12. These rules and strict-root tolerances are fixed.

Selection uses only controlled terminations with no numerical error, finite
coefficients, finite original-matrix training predictors/risk/loss, and relative
coordinate reconstruction error at most 1e-8. The minimum recomputed training
loss wins. Loss differences within `1e-12 * max(1, abs(minimum_loss))` use priority
cell mean, global mean, perturbed cell. No new data or held-out outcomes are
accepted by the fitting API or used for selection. Prediction failures never
trigger reselection.

`fit$selection$table` retains eligibility, loss, termination and selected flags
for every start. `fit$endpoints` retains each coefficient vector, diagnostic
metrics, accepted path, proposals, trace and error. `fit$selection$spread`
measures the eligible training loss range; it does not establish prediction
stability. Users can inspect predictions from retained coefficients to assess
start sensitivity on a prespecified support.

A selected endpoint is reported as `FINITE_ENDPOINT_SELECTED`; if all endpoints
fail eligibility, the fit reports `NO_ELIGIBLE_ENDPOINT` and prediction errors.
Numerical failures are retained as `NUMERIC_ERROR`, with error text. Severe
column conditioning can make numerical solves fail even when all inputs are
finite. Such failures are explicit and do not select a substitute method.

## Prediction and scientific limits

`predict(fit, newdata, type = "risk")` also accepts `"survival"` and `"link"`.
Named new-data columns are matched and reordered to training names; unnamed
matrices use positional columns. Output follows new-data row order and names.
Missing predictor rows return NA. Infinite predictors are errors; overflow of
a finite matrix product returns NA for the affected rows. Finite predictors can
produce exact probability zero or one through floating-point saturation.
Zero-row new data work. The caller must preserve cell encoding and fitted
preprocessing; dimensions and names are validated, while extrapolation and
new-data cell geometry are not screened.

A finite stopped endpoint need not solve the estimating equation or attain a
local or global minimum. Strict-root, derivative rank, full-update and Hessian
diagnostics are retained but do not gate finite predictions. A successful
synthetic test is not validation of clinical performance, population accuracy
or inferential validity.

This package is not a general GEE replacement and is not Cox partial likelihood.
It supplies no arbitrary weights, families, working correlations, penalties,
monotonicity repair, robust standard errors, confidence intervals or coefficient
inference. The optimization Hessian is not a robust covariance estimator.

## Validation and source provenance

`tests/synthetic.R` runs with base R against the installed package and uses only
synthetic inputs. It checks derivative finite differences, trust-step conditions,
known-mean recovery, raw and boundary responses, deterministic starts, selection
rules, finite unconverged endpoints, explicit numerical failures, design/input
rejection, and prediction naming, missingness and invariance of selection.
Runnable synthetic examples are included in the help pages. Local reference
parity checks and check logs are held outside the source distribution.

The generic mathematics, trust solver, starts and selection were transcribed
from the existing independently developed project implementation. No geepack
source is included or called, and no private data, result vectors or research
paths are needed. The public package contains code, documentation and synthetic
tests only; research datasets and clinical results are not distributed.
