common_starts <- function(d) {
  ix<-d$intercepts;stopifnot(all(ix%in%seq_len(ncol(d$x))),all(rowSums(d$x[,ix,drop=FALSE])==1))
  clip<-function(v)pmin(1-1e-6,pmax(1e-6,v)); link<-function(v)log(-log1p(-clip(v)))
  b1<-b2<-rep(0,ncol(d$x));b1[ix]<-link(mean(d$y))
  for(j in ix)b2[j]<-link(mean(d$y[d$x[,j]==1]))
  stopifnot(all(is.finite(c(b1,b2))))
  list(global_mean=b1,cell_mean=b2)
}
