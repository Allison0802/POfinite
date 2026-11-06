# Pure unchanged-objective mathematics, no data or package dependencies.
link_values <- function(eta) {
  a <- exp(eta); mu <- -expm1(-a)
  d <- exp(eta-a); d2 <- numeric(length(eta)); ok <- is.finite(a)
  d2[ok] <- d[ok]*(1-a[ok])
  list(mu=mu,d=d,d2=d2)
}
orthogonalize <- function(x) {
  sx<-sqrt(colMeans(x*x)); q<-qr(sweep(x,2,sx,'/'),tol=1e-12)
  stopifnot(q$rank==ncol(x)); tmp<-matrix(0,ncol(x),ncol(x))
  tmp[q$pivot,]<-backsolve(qr.R(q),diag(ncol(x)))
  T<-sweep(tmp,1,sx,'/')*sqrt(nrow(x)); z<-qr.Q(q)*sqrt(nrow(x))
  stopifnot(max(abs(x%*%T-z))<1e-8)
  list(T=T,z=z)
}
state <- function(z,y,theta,hessian=FALSE) {
  eta<-drop(z%*%theta); v<-link_values(eta); r<-y-v$mu; n<-nrow(z)
  g<-drop(crossprod(z,-v$d*r))/n
  H<-if(hessian)crossprod(z,z*(v$d^2-r*v$d2))/n else NULL
  c(v,list(eta=eta,r=r,L=sum(r*r)/(2*n),g=g,H=H))
}
tr_step <- function(g,H,radius) {
  e<-eigen(H,symmetric=TRUE); q<-drop(crossprod(e$vectors,g)); h<-e$values
  tol<-max(1e-14,max(abs(h))*1e-12); lo<-max(0,-min(h))
  # Regular positive definite Newton point inside the ball.
  if(min(h)>tol){p<- -drop(e$vectors%*%(q/h));if(sqrt(sum(p*p))<=radius)return(list(p=p,lambda=0,hard=FALSE))}
  den<-h+lo; use<-den>tol; w<-numeric(length(g)); w[use]<- -q[use]/den[use]
  if(all(abs(q[!use])<=tol)&&sqrt(sum(w*w))<=radius){
    # Hard case, or flat positive-semidefinite subproblem. Choose sign giving descent.
    if(lo>0){j<-which.min(h);w[j]<-if(q[j]>0)-sqrt(max(0,radius^2-sum(w*w)))else sqrt(max(0,radius^2-sum(w*w)))}
    return(list(p=drop(e$vectors%*%w),lambda=lo,hard=TRUE))
  }
  lower<-lo; upper<-max(1,lo+1)
  norm_at<-function(lam)sqrt(sum((q/(h+lam))^2))
  while(norm_at(upper)>radius)upper<-2*upper+1
  for(i in 1:100){mid<-(lower+upper)/2;if(norm_at(mid)>radius)lower<-mid else upper<-mid}
  list(p= -drop(e$vectors%*%(q/(h+upper))),lambda=upper,hard=FALSE)
}
endpoint <- function(x,y,z,T,b) {
  n<-nrow(x); eta<-drop(x%*%b); v<-link_values(eta); r<-y-v$mu
  J<-z*v$d; ss<-svd(J,nu=0); rank<-sum(ss$d>max(ss$d)*1e-10)
  H<-crossprod(z,z*(v$d^2-r*v$d2))/n; ev<-eigen(H,symmetric=TRUE,only.values=TRUE)$values
  gn<-rep(NA_real_,ncol(x)); pred<-NA_real_
  if(rank==ncol(x)){
    # QR full update independent of the SVD rank diagnostic.
    gn<-drop(T%*%qr.coef(qr(J,tol=1e-12),r))
    pred<-max(abs(link_values(drop(x%*%(b+gn)))$mu-v$mu))
  }
  score<-max(abs(crossprod(sweep(x,2,sqrt(colMeans(x*x)),'/'),v$d*r)))/n
  theta<-solve(T,b); rec<-drop(z%*%theta)
  vals<-c(L=sum(r*r)/(2*n),mse=mean(r*r),scaled_score=score,p=ncol(x),rank=rank,
    sigma_min=min(ss$d),sigma_next=sort(ss$d)[2],sigma_max=max(ss$d),condition=max(ss$d)/min(ss$d),
    hessian_min=min(ev),hessian_max=max(ev),negative_eigen=sum(ev< -max(abs(ev))*1e-8),
    full_step_max=max(abs(gn)),full_step_prediction_change=pred,
    eta_min=min(eta),eta_max=max(eta),derivative_min=min(v$d),derivative_max=max(v$d),
    beta_l2=sqrt(sum(b*b)),beta_max=max(abs(b)),
    coordinate_abs_error=max(abs(eta-rec)),coordinate_relative_error=max(abs(eta-rec)/pmax(1,abs(eta))),
    strict_pass=rank==ncol(x)&&score<=1e-10&&is.finite(max(abs(gn)))&&max(abs(gn))<=1e-8&&is.finite(pred)&&pred<=1e-9&&min(ev)>=-max(abs(ev))*1e-8)
  list(values=vals,sv=ss,eta=eta,mu=v$mu)
}
