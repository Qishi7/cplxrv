library(cplxshrink)
X<-matrix(c(-1,0,1,2),ncol=1);y<-c(-.5,.1,.7,1)+1i*c(.2,-.1,.1,-.2)
for(f in c("student","lasso","gdp"))for(d in c(FALSE,TRUE)) {
  set.seed(25);a<-cplxshrink:::smcg_sample(if(d)y else Re(y),X,family=f,complex_model=d,
                             rho=if(d).2+.3i else 0,iterations=30,warm=3,thin=2)
  set.seed(25);b<-cplxshrink:::smcg_sample(if(d)y else Re(y),X,family=f,complex_model=d,
                             rho=if(d).2+.3i else 0,iterations=30,warm=3,thin=2)
  stopifnot(identical(a,b),all(is.finite(a$draws)),identical(a$sampleidx,seq(5L,29L,2L)))
}
