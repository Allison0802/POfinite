# Matrix validation is explicit: preprocessing and aliases belong to the caller.
.pf_matrix <- function(x, what) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || ncol(x) < 1L)
    stop(what, " must be a real numeric matrix with at least one column")
  nm <- colnames(x)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
    stop(what, " column names must be nonempty, unique and nonmissing")
  x
}

pofinite_fit <- function(x, y, intercepts) {
  x <- .pf_matrix(x, "x")
  if (!nrow(x) || any(!is.finite(x))) stop("x must have finite training rows")
  if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) ||
      length(y) != nrow(x) || any(!is.finite(y)))
    stop("y must be a finite numeric vector with one value per training row")
  if (!is.numeric(intercepts) || !length(intercepts) || anyNA(intercepts) ||
      !identical(as.numeric(intercepts), as.numeric(seq_along(intercepts))))
    stop("intercepts must be the consecutive leading column indices 1:k")
  k <- length(intercepts)
  if (k >= ncol(x)) stop("unsupported design: the frozen perturbation requires p > k")
  ix <- seq_len(k)
  cells <- x[, ix, drop = FALSE]
  if (any(!cells %in% c(0, 1)) || any(rowSums(cells) != 1) || any(colSums(cells) == 0))
    stop("intercept columns must be nonempty disjoint 0/1 cells partitioning the rows")
  sx <- sqrt(colMeans(x * x))
  if (nrow(x) < ncol(x) || any(!is.finite(sx)) || any(sx == 0))
    stop("unsupported design: finite positive column scales and n >= p required")
  q <- qr(sweep(x, 2, sx, "/"), tol = 1e-12)
  if (q$rank != ncol(x)) stop("unsupported design: remove exact aliases before fitting")
  if (!identical(q$pivot[ix], ix)) stop("unsupported design: QR must preserve leading intercepts")
  d <- list(x = x, y = y, intercepts = ix)
  starts <- tryCatch(common_starts(d), error = function(e)
    stop("unsupported design: frozen start geometry unavailable (",
         conditionMessage(e), ")", call. = FALSE))
  fits <- lapply(starts, function(b) {
    tryCatch(common_solve_fast(d, b), error = function(e)
      list(beta = rep(NA_real_, ncol(x)), metrics = NULL, reason = "NUMERIC_ERROR",
           error = conditionMessage(e), accepted = NA_integer_, trace = data.frame(),
           proposals = list(), beta_path = list(b)))
  })
  selection <- select_training(d, fits)
  selected <- selection$selected
  selection$table$selected <- !is.na(selected) & selection$table$start == selected
  ok <- !is.na(selected)
  b <- if (ok) fits[[selected]]$beta else rep(NA_real_, ncol(x))
  names(b) <- colnames(x)
  structure(list(coefficients = b, selected = selected, selection = selection,
                 endpoints = fits, starts = starts, intercepts = ix,
                 column_names = colnames(x), n_columns = ncol(x), n = nrow(x),
                 status = if (ok) "FINITE_ENDPOINT_SELECTED" else "NO_ELIGIBLE_ENDPOINT"),
            class = "pofinite")
}

predict.pofinite <- function(object, newdata, type = c("risk", "survival", "link"), ...) {
  type <- match.arg(type)
  if (!identical(object$status, "FINITE_ENDPOINT_SELECTED"))
    stop("prediction unavailable: no eligible training endpoint")
  x <- .pf_matrix(newdata, "newdata")
  if (ncol(x) != object$n_columns) stop("newdata has the wrong number of columns")
  expected <- object$column_names
  actual <- colnames(x)
  if (is.null(expected) != is.null(actual)) stop("newdata naming must match the training matrix")
  if (!is.null(expected)) {
    if (!setequal(expected, actual)) stop("newdata column names differ from training")
    x <- x[, expected, drop = FALSE]
  }
  if (any(is.infinite(x))) stop("newdata cannot contain infinite predictors")
  answer <- rep(NA_real_, nrow(x))
  valid <- rowSums(is.na(x)) == 0L
  eta <- drop(x[valid, , drop = FALSE] %*% object$coefficients)
  good <- is.finite(eta)
  values <- rep(NA_real_, length(eta))
  values[good] <- switch(type, link = eta[good], risk = -expm1(-exp(eta[good])),
                         survival = exp(-exp(eta[good])))
  answer[valid] <- values
  names(answer) <- rownames(x)
  answer
}

coef.pofinite <- function(object, ...) object$coefficients

print.pofinite <- function(x, ...) {
  cat("Finite-iteration pseudo-observation predictor\n")
  cat("Status:", x$status, "\n")
  cat("Training rows:", x$n, "; columns:", x$n_columns, "\n")
  if (!is.na(x$selected)) cat("Selected training start:", x$selected, "\n")
  print(x$selection$table, row.names = FALSE)
  invisible(x)
}
