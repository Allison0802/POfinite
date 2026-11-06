# Finite-iteration exact-Hessian trust-region solver.
common_solve_fast <- function(d,b0,checkpoint=function(value)NULL,resume=NULL) {
  x<-d$x;y<-d$y;ot<-orthogonalize(x);z<-ot$z;T<-ot$T
  theta<-solve(T,b0);radius<-.5;accepted<-0L;bs<-list(b0);rows<-list();proposals<-list();reason<-'TRUST_PROPOSAL_LIMIT'
  if(!is.null(resume)){theta<-resume$theta;radius<-resume$radius;accepted<-resume$accepted;bs<-resume$beta_path;rows<-resume$trace;proposals<-resume$proposals}
  first<-length(rows)+1L
  if(accepted>=400L)reason<-'TRUST_ACCEPTED_LIMIT'
  if(radius<1e-12)reason<-'TRUST_RADIUS_STALL'
  error_message<-NULL;xs<-sweep(x,2,sqrt(colMeans(x*x)),'/')
  tryCatch(if(first<=600L&&accepted<400L&&radius>=1e-12)for(it in seq.int(first,600L)) {
    if(any(!is.finite(z%*%theta))||any(!is.finite(x%*%(T%*%theta))))stop('NONFINITE_TRAIN_ETA')
    v<-state(z,y,theta,TRUE)
    ov<-link_values(drop(x%*%drop(T%*%theta)))
    exact_score<-max(abs(crossprod(xs,ov$d*(y-ov$mu))))/nrow(x)
    if(exact_score<=1e-10){mm<-endpoint(x,y,z,T,drop(T%*%theta))$values;if(mm[['strict_pass']]==1){reason<-'STRICT_ROOT';break}}
    s<-tr_step(v$g,v$H,radius);step<-s$p;w<-state(z,y,theta+step)
    predicted<- -sum(v$g*step)-sum(step*drop(v$H%*%step))/2
    reduction<- -mean((w$mu-v$mu)*(.5*(w$mu+v$mu)-y))
    rho<-if(predicted>0)reduction/predicted else -Inf
    take<-is.finite(rho)&&rho>.1&&reduction>0
    rows[[it]]<-data.frame(iter=it,L_before=v$L,L_after=w$L,accepted=take,reduction=reduction,predicted=predicted,rho=rho,radius=radius,lambda=s$lambda)
    proposals[[it]]<-list(theta=theta,step=step)
    pn<-sqrt(sum(step*step));if(!is.finite(rho)||rho<.25)radius<-radius/4 else if(rho>.75&&pn>=.9*radius)radius<-min(8,2*radius)
    if(take){theta<-theta+step;accepted<-accepted+1L;bs[[length(bs)+1]]<-drop(T%*%theta)}
    if(it%%25L==0L)checkpoint(list(theta=theta,radius=radius,accepted=accepted,beta_path=bs,trace=rows,proposals=proposals,iteration=it))
    if(radius<1e-12){reason<-'TRUST_RADIUS_STALL';break}
    if(accepted>=400){reason<-'TRUST_ACCEPTED_LIMIT';break}
  },error=function(e){reason<<-'NUMERIC_ERROR';error_message<<-conditionMessage(e)})
  b<-drop(T%*%theta);mm<-endpoint(x,y,z,T,b)$values;if(mm[['strict_pass']]==1)reason<-'STRICT_ROOT'
  list(beta=b,metrics=mm,reason=reason,error=error_message,accepted=accepted,trace=if(length(rows))do.call(rbind,rows)else data.frame(),proposals=proposals,beta_path=bs)
}
