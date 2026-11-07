test_that("dimensions, proper priors, schedules and fixed relation are enforced", {
  for (p in c(1L,3L)) for (model in c("real","complex")) for (prior in c("student","lasso","gdp")) {
    X <- matrix(seq(-1,1,length.out=8*p),8,p)
    y <- seq(-.4,.5,length.out=8); if(model=="complex") y<-y+1i*rev(y)
    z <- cplx_shrink(y,X,model=model,prior=prior,iterations=20,warm=3,thin=2)
    expect_equal(dim(z$beta),c(8,if(model=="complex")2*p else p))
    expect_equal(z$sampleidx,seq(5,19,2))
    expect_true(all(is.finite(z$draws)))
    expect_true(all(c(z$w,z$s,z$hyper)>0))
  }
  X<-matrix(1:8,8,1); y<-seq_len(8)
  expect_error(cplx_shrink(y,X,relation=1),"rho")
  expect_error(cplx_shrink(y,X,model="real",relation=.1),"rho")
  expect_error(cplx_shrink(y+1i,X,model="real"),"real data")
  expect_error(cplx_shrink(y[-1],X),"data")
  expect_error(cplx_shrink(y,X,warm=10,iterations=10),"iterations")
  expect_error(cplx_shrink(y,X,hyper=list(nu=-1)),"positive")
  expect_error(cplx_shrink(y,X,hyper=list(learn_rho=TRUE)),"hyper")
  expect_error(cplx_shrink(y,X,start=list(L=1)),"unused")
  expect_error(cplxshrink:::smcg_sample(y,X,rho_b=.1),"rho_b")
  expect_error(cplxshrink:::smcg_sample(y,X,learn_rho=TRUE),"disabled")
  expect_error(cplxshrink:::smcg_sample(y,X,griddy=TRUE),"disabled")
  expect_equal(cplxshrink:::covariance_cpp(2L,2L,0,0)$C,diag(4))
})

test_that("complex and structured-real likelihood and posterior parameters agree", {
  X<-matrix(c(1,2,3,4,-1,1),3,2)+1i*matrix(c(.1,-.2,.4,.5,.3,-.1),3,2)
  b<-c(.3+.2i,-.1+.4i);y<-c(.2+.3i,.4-.2i,-.5+.1i);rho<-.3-.2i;s<-.7
  H<-rbind(cbind(Re(X),-Im(X)),cbind(Im(X),Re(X)))
  yy<-c(Re(y),Im(y));bb<-c(Re(b),Im(b))
  expect_equal(drop(H%*%bb),c(Re(X%*%b),Im(X%*%b)),tolerance=1e-14)
  e<-y-drop(X%*%b);ec<-sum((Mod(e)^2-Re(rho*Conj(e)^2))/(1-Mod(rho)^2))/s
  K<-cplxshrink:::covariance_cpp(3L,2L,Re(rho),Im(rho))$K
  expect_equal(drop(crossprod(yy-H%*%bb,K%*%(yy-H%*%bb)))/s,ec,tolerance=1e-13)
  w<-c(.5,2);A<-crossprod(H,K%*%H)+diag(rep(1/w,2))
  h<-crossprod(H,K%*%yy)
  # Direct independent whitening produces the same analytic posterior.
  C<-solve(K);L<-t(chol(C));HW<-forwardsolve(L,H);yw<-forwardsolve(L,yy)
  Aw<-crossprod(HW)+diag(rep(1/w,2))
  expect_equal(solve(A,h),solve(Aw,crossprod(HW,yw)),tolerance=1e-13)
  expect_equal(s*solve(A),s*solve(Aw),tolerance=1e-13)
  for(f in c("student","lasso","gdp")) {
    set.seed(142);z<-cplx_shrink(y,X,prior=f,relation=rho,iterations=100,warm=10,start=list(w=w,s=s))
    set.seed(142);q<-cplxshrink:::canonical_cpp(yy,H,2L,f,Re(rho),Im(rho),100L,10L,1L,rep(0,4),w,s,1,rep(1,2),1,1,.5,.5,1,1,.5,.5)
    expect_identical(unname(z$beta),q$beta)
    expect_identical(z$s,q$s)
  }
})

test_that("lasso global conditional uses dimension-consistent gamma shape", {
  for(d in 1:2) for(p in c(1L,3L)) {
    X<-matrix(seq(-1,1,length.out=8*p),8,p);y<-seq(-.4,.5,length.out=8)
    if(d==2)y<-y+1i*rev(y)
    set.seed(27000+10*d+p)
    z<-cplx_shrink(y,X,prior="lasso",model=if(d==2)"complex" else "real",
      hyper=list(r=1.3,delta=.7),iterations=6000,warm=1000)
    U<-pgamma(as.vector(z$hyper),shape=(d+1)*p/2+1.3,rate=.7+rowSums(z$w)/2)
    expect_lt(unname(stats::ks.test(U,"punif")$statistic),.035)
  }
})

test_that("GDP joint restoration has the correct collapsed gamma conditional", {
  for(d in 1:2) {
    set.seed(5600+d);z<-cplxshrink:::gdp_restore_cpp(20000L,.8,.7,6,1.2,d)
    u<-pgamma(z[,1],shape=6+d,rate=1.2+sqrt(.8/.7))
    expect_lt(unname(stats::ks.test(u,"punif")$statistic),.018)
    expect_true(all(is.finite(z)&z>0))
  }
  # At chi=0 the local conditional is Gamma(1/2, psi/2), not an atom.
  set.seed(778);w<-cplxshrink:::local_w_cpp(20000L,0,2)
  expect_equal(mean(w),.5,tolerance=.025)
})
