base_starts<-common_starts
common_starts<-function(d){
 ss<-base_starts(d);x<-d$x;k<-length(d$intercepts);p<-ncol(x);n<-nrow(x)
 stopifnot(identical(as.integer(d$intercepts),seq_len(k)),p>k)
 q<-qr(sweep(x,2,sqrt(colMeans(x*x)),'/'),tol=1e-12);stopifnot(identical(q$pivot[seq_len(k)],seq_len(k)))
 ot<-orthogonalize(x);v<-c(rep(0,k),rep(c(1,-1),length.out=p-k)*.25/sqrt(p-k));inc<-drop(ot$T%*%v);eta<-drop(x%*%inc)
 stopifnot(abs(sqrt(mean(eta^2))-.25)<1e-10,max(abs(crossprod(x[,d$intercepts,drop=FALSE],eta)/colSums(x[,d$intercepts,drop=FALSE])))<1e-10,max(abs(eta-drop(ot$z%*%v)))<1e-8)
 ss$perturbed_cell<-ss$cell_mean+inc;ss
}
select_training<-function(d,fits,omit=character()){
 pri<-c('cell_mean','global_mean','perturbed_cell');rows<-list();ot<-orthogonalize(d$x)
 for(n in names(fits)){
  f<-fits[[n]];b<-f$beta;finite<-length(b)==ncol(d$x)&&all(is.finite(b));L<-coord<-NA_real_
  if(finite){et<-drop(d$x%*%b);mu<- -expm1(-exp(et));finite<-all(is.finite(c(et,mu)));if(finite){L<-mean((d$y-mu)^2)/2;coord<-max(abs(et-drop(ot$z%*%solve(ot$T,b)))/pmax(1,abs(et)))}}
  eligible<-is.null(f$error)&&finite&&is.finite(L)&&is.finite(coord)&&coord<=1e-8&&f$reason%in%c('STRICT_ROOT','TRUST_RADIUS_STALL','TRUST_PROPOSAL_LIMIT','TRUST_ACCEPTED_LIMIT')&&!(n%in%omit)
  rows[[n]]<-data.frame(start=n,loss=L,training_coordinate_error=coord,eligible=eligible,reason=f$reason)
 };tab<-do.call(rbind,rows);at<-which(tab$eligible);chosen<-NA_character_
 if(length(at)){lo<-min(tab$loss[at]);ties<-tab$start[at][tab$loss[at]<=lo+1e-12*max(1,abs(lo))];chosen<-pri[pri%in%ties][1]}
 list(selected=chosen,table=tab,spread=if(length(at))diff(range(tab$loss[at]))else NA_real_)
}
