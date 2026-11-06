library(pofinite)

# Self-contained synthetic tests: no private data, reference files or network.
nchecks <- 0L
check <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  nchecks <<- nchecks + 1L
  cat(sprintf("ok %02d - %s\n", nchecks, label))
}
near <- function(a,b,tol=1e-9) length(a)==length(b) && all(is.finite(c(a,b))) && max(abs(a-b),0)<tol
fails <- function(expr, text) {
  e <- tryCatch({force(expr); NULL}, error=identity)
  inherits(e,"error") && grepl(text,conditionMessage(e),fixed=TRUE)
}
internal <- function(name) getFromNamespace(name,"pofinite")
link <- function(eta) -expm1(-exp(eta))
x <- cbind(cellA=rep(c(1,0),each=12), cellB=rep(c(0,1),each=12),
           z=rep(seq(-1,1,length.out=12),2))
y <- c(seq(.05,.7,length.out=12),seq(.2,.9,length.out=12))
fit <- pofinite_fit(x,y,1:2)
check("fit has selected and nonselected endpoint flags", inherits(fit,"pofinite") &&
      sum(fit$selection$table$selected)==1 && nrow(fit$selection$table)==3)
check("synthetic loss regression and tie priority", fit$selected=="cell_mean" &&
      near(fit$selection$table$loss, rep(0.0003209104,3),1e-10))
check("training original-matrix loss is recomputed", near(
      fit$selection$table$loss, vapply(fit$endpoints,function(f)
        mean((y-link(drop(x%*%f$beta)))^2)/2,numeric(1)),1e-16))
check("fit and starts are deterministic", identical(fit,pofinite_fit(x,y,1:2)))
check("coefficient extraction retains training names", identical(coef(fit),fit$coefficients) &&
      identical(names(coef(fit)),colnames(x)))
check("print returns the fit invisibly", {
  visible <- NULL
  capture.output(visible <- withVisible(print(fit)))
  identical(visible$value,fit) && !visible$visible
})

# Derivatives are checked independently with central finite differences.
ot <- internal("orthogonalize")(x)
state <- internal("state")
theta <- c(-.4,.2,.1)
h <- 1e-5
v <- state(ot$z,y,theta,TRUE)
g <- vapply(seq_along(theta),function(j) {
  delta <- replace(rep(0,length(theta)),j,h)
  (state(ot$z,y,theta+delta)$L-state(ot$z,y,theta-delta)$L)/(2*h)
},numeric(1))
H <- vapply(seq_along(theta),function(j) {
  delta <- replace(rep(0,length(theta)),j,h)
  (state(ot$z,y,theta+delta)$g-state(ot$z,y,theta-delta)$g)/(2*h)
},numeric(length(theta)))
check("analytic gradient matches finite differences",near(g,v$g,1e-9))
check("exact Hessian matches gradient differences",near(H,v$H,1e-9))
check("orthogonal training coordinates reconstruct",near(x%*%ot$T,ot$z,1e-12) &&
      near(crossprod(ot$z)/nrow(x),diag(ncol(x)),1e-12))
tr <- internal("tr_step")
s <- tr(c(.1,.2),diag(c(2,3)),1)
check("interior positive-definite trust solution",near(s$p,c(-.05,-.2/3)) && s$lambda==0)
s <- tr(c(.1,.2),diag(c(-2,3)),.5)
check("indefinite trust step satisfies KKT and radius",near(sqrt(sum(s$p^2)),.5) &&
      near(drop((diag(c(-2,3))+diag(s$lambda,2))%*%s$p),-c(.1,.2)) && s$lambda>=2)
s <- tr(c(0,0),diag(c(-2,3)),.5)
check("hard-case trust step follows negative curvature",near(sum(s$p^2),.25) &&
      s$hard && near(s$lambda,2))

# Start clipping affects only initialization, never the raw response objective.
inc <- fit$starts$perturbed_cell-fit$starts$cell_mean
check("perturbation has prescribed RMS and zero cell means",near(sqrt(mean((x%*%inc)^2)),.25) &&
      max(abs(crossprod(x[,1:2],x%*%inc)))<1e-12)
x1 <- cbind(cell=1,z=seq(-1,1,length.out=24))
outside <- rep(c(-.4,0,.3,.8,1,1.4),4)
fout <- pofinite_fit(x1,outside,1L)
check("outside-range responses retain raw loss",near(
      fout$selection$table$loss[fout$selection$table$selected],
      mean((outside-predict(fout,x1))^2)/2,1e-14))
known <- c(-.7,.5)
fknown <- pofinite_fit(x1,link(drop(x1%*%known)),1L)
check("known noiseless mean and coefficients recover",near(coef(fknown),known,1e-7) &&
      near(predict(fknown,x1),link(drop(x1%*%known)),1e-8))
fmean <- pofinite_fit(x1,rep(.3,nrow(x1)),1L)
check("known constant mean recovers with extra predictor",near(predict(fmean,x1),rep(.3,nrow(x1)),1e-9))
fzero <- pofinite_fit(x1,rep(0,nrow(x1)),1L)
fone <- pofinite_fit(x1,rep(1,nrow(x1)),1L)
check("boundary raw responses give finite probabilities",all(is.finite(c(predict(fzero,x1),predict(fone,x1)))))
check("unconverged finite endpoints remain eligible",all(fzero$selection$table$eligible) &&
      all(vapply(fzero$endpoints,function(f)f$metrics[['strict_pass']]==0,logical(1))) &&
      !is.na(fzero$selected))
check("all controlled termination categories alone do not gate selection", {
  d <- list(x=x,y=y)
  fs <- fit$endpoints
  fs[[1]]$reason <- 'TRUST_PROPOSAL_LIMIT'
  fs[[2]]$reason <- 'TRUST_ACCEPTED_LIMIT'
  fs[[3]]$reason <- 'TRUST_RADIUS_STALL'
  all(internal('select_training')(d,fs)$table$eligible)
})
check("numerical-error endpoints are excluded despite finite coefficients", {
  fs <- fit$endpoints
  fs$cell_mean$error <- 'synthetic numerical error'
  fs$cell_mean$reason <- 'NUMERIC_ERROR'
  sel <- internal('select_training')(list(x=x,y=y),fs)
  !sel$table$eligible[sel$table$start=='cell_mean'] && sel$selected!='cell_mean'
})
check("unknown termination is excluded", {
  fs <- fit$endpoints
  fs$cell_mean$reason <- 'UNKNOWN'
  !internal('select_training')(list(x=x,y=y),fs)$table$eligible[2]
})
check("all numeric errors yield explicit failure status", {
  # Fault injection exercises wrapper error reporting without changing solver code.
  f <- pofinite_fit
  environment(f) <- list2env(list(common_solve_fast=function(...)stop('synthetic fault')),
                             parent=asNamespace('pofinite'))
  failed <- f(x,y,1:2)
  failed$status=='NO_ELIGIBLE_ENDPOINT' && is.na(failed$selected) &&
    all(is.na(coef(failed))) && all(vapply(failed$endpoints,function(e)
      identical(e$reason,'NUMERIC_ERROR') && identical(e$error,'synthetic fault'),logical(1))) &&
    fails(predict(failed,x),'no eligible training endpoint')
})
check("training loss sensitivity and omitted-start selection retained", {
  fs <- fit$endpoints
  fs$perturbed_cell$beta <- fit$starts$perturbed_cell
  sel <- internal('select_training')(list(x=x,y=y),fs)
  loo <- internal('select_training')(list(x=x,y=y),fs,omit='cell_mean')
  sel$spread>1e-3 && sel$selected=='cell_mean' && loo$selected=='global_mean'
})

# Prediction uses frozen coefficients, with explicit row and column semantics.
newx <- x[c(3,1,2),c(3,2,1)]
rownames(newx) <- c('third','first','second')
p <- predict(fit,newx)
check("named predictor order and output row names are preserved",identical(names(p),rownames(newx)) &&
      near(unname(p),predict(fit,x[c(3,1,2),])))
check("risk survival and link predictions agree",near(predict(fit,x)+predict(fit,x,type='survival'),rep(1,nrow(x))) &&
      near(predict(fit,x),link(predict(fit,x,type='link'))))
newx[2,1] <- NA_real_
newx[3,2] <- NaN
check("NA and NaN predictor rows remain in place",identical(unname(is.na(predict(fit,newx))),c(FALSE,TRUE,TRUE)) &&
      identical(names(predict(fit,newx)),rownames(newx)))
check("zero-row prediction supported",length(predict(fit,x[FALSE,,drop=FALSE]))==0L)
check("prediction never changes selection or endpoints", {
  before <- serialize(fit,NULL)
  invisible(predict(fit,x*2))
  identical(before,serialize(fit,NULL)) &&
    !any(c('newdata','test_y','holdout') %in% names(formals(pofinite_fit)))
})
unnamed <- unname(x1)
fu <- pofinite_fit(unnamed,rep(.3,nrow(unnamed)),1L)
check("unnamed matrix uses positional columns",near(predict(fu,unnamed),rep(.3,nrow(unnamed))))
check("nonfinite prediction products return missing rows", {
  ff <- fit
  ff$coefficients[] <- 2
  xx <- matrix(1e308,1,3,dimnames=list(NULL,colnames(x)))
  is.na(predict(ff,xx))
})

# Unsupported inputs fail explicitly instead of silently changing the estimator.
check("intercept-only design rejected",fails(pofinite_fit(matrix(1,5,1),rep(.2,5),1L),'p > k'))
check("exact aliases rejected",fails(pofinite_fit(cbind(a=1,b=1,z=1:5),rep(.2,5),1L),'aliases'))
check("training missingness and response shape rejected",fails(pofinite_fit(x,replace(y,1,NA),1:2),'finite numeric') &&
      fails(pofinite_fit(x,matrix(y),1:2),'finite numeric') &&
      fails(pofinite_fit(replace(x,1,NA),y,1:2),'finite training'))
check("bad intercept positions and cell geometry rejected",fails(pofinite_fit(x,y,2:3),'leading') &&
      fails(pofinite_fit(replace(x,1,.5),y,1:2),'0/1 cells'))
check("nonmatrix and complex inputs rejected",fails(pofinite_fit(as.data.frame(x),y,1:2),'real numeric matrix') &&
      fails(pofinite_fit(x+0i,y,1:2),'real numeric matrix'))
check("zero scale and underdetermined matrices rejected",fails(pofinite_fit(cbind(cell=rep(1,5),z=rep(0,5)),rep(.2,5),1L),'column scales') &&
      fails(pofinite_fit(x[c(1,13),],y[c(1,13)],1:2),'n >= p'))
check("duplicate or missing training column names rejected", {
  xx <- x; colnames(xx)[2] <- colnames(xx)[1]
  xx2 <- x; colnames(xx2)[2] <- NA_character_
  fails(pofinite_fit(xx,y,1:2),'column names') && fails(pofinite_fit(xx2,y,1:2),'column names')
})
check("new-data infinite, wrong dimension and name mismatch rejected", {
  xx <- x; xx[1,1] <- Inf
  xx2 <- x; colnames(xx2)[1] <- 'other'
  fails(predict(fit,xx),'infinite') && fails(predict(fit,x[,1:2]),'number of columns') &&
    fails(predict(fit,xx2),'column names differ') && fails(predict(fit,unname(x)),'naming') &&
    fails(predict(fu,x1),'naming')
})
check("unsupported API extensions fail",fails(pofinite_fit(x,y,1:2,weights=rep(1,nrow(x))),'unused argument') &&
      fails(predict(fit,x,type='response'),'arg'))
check("severe column scaling gives explicit endpoint failures", {
  results <- lapply(c(1e16,1e20,1e-20),function(scale) {
    xx <- x1; xx[,2] <- xx[,2]*scale
    f <- pofinite_fit(xx,rep(.3,nrow(xx)),1L)
    f$status=='NO_ELIGIBLE_ENDPOINT' && all(vapply(f$endpoints,function(e)
      identical(e$reason,'NUMERIC_ERROR') && !is.null(e$error),logical(1)))
  })
  all(unlist(results))
})
check("near-collinear retained design reports finite or explicit failure status", {
  xx <- cbind(x1,near=x1[,2]+1e-8*sin(seq_len(nrow(x1))))
  f <- tryCatch(pofinite_fit(xx,link(drop(x1%*%c(-.7,.5))),1L),error=identity)
  if (inherits(f,'error')) grepl('unsupported design:',conditionMessage(f),fixed=TRUE) else
    f$status %in% c('FINITE_ENDPOINT_SELECTED','NO_ELIGIBLE_ENDPOINT') &&
      (f$status!='FINITE_ENDPOINT_SELECTED' || all(is.finite(predict(f,xx))))
})
check("extreme finite responses fail explicitly when loss overflows", {
  f <- suppressWarnings(pofinite_fit(x1,rep(1e308,nrow(x1)),1L))
  f$status=='NO_ELIGIBLE_ENDPOINT' && !any(f$selection$table$eligible)
})
cat(sprintf('PASS: %d synthetic checks\n',nchecks))
