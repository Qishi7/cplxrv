test_that("six samplers reproduce independent posterior quadrature moments", {
  skip_if_not_installed("jsonlite")
  data<-jsonlite::fromJSON(test_path("fixtures/reference_data.json"))
  ref<-jsonlite::fromJSON(test_path("fixtures/reference_quadrature.json"),simplifyVector=FALSE)
  for(model in c("real","complex"))for(prior in c("student","lasso","gdp")) {
    yy<-if(model=="real")data$yr else data$yr+1i*data$yi
    hyp<-as.list(data$parameters)
    use<-c("a_s","b_s",switch(prior,student=if(model=="complex")"nu" else c("a_tau","b_tau"),lasso=c("r","delta"),gdp=c("alpha","eta")))
    chains<-lapply(1:4,function(k){
      set.seed(252500+100*match(paste(model,prior),c("real student","real lasso","real gdp","complex student","complex lasso","complex gdp"))+k)
      z<-cplx_shrink(yy,matrix(data$x,ncol=1),prior=prior,model=model,
        relation=if(model=="complex")data$rho_re+1i*data$rho_im else 0,
        hyper=hyp[use],iterations=8000,warm=2000,start=list(s=c(.1,.4,1,2)[k],w=c(.2,1,3,5)[k]))
      mm<-cbind(z$beta,s=as.vector(z$s),w=as.vector(z$w))
      colnames(mm)[seq_len(ncol(z$beta))]<-paste0("b",seq_len(ncol(z$beta)))
      for(j in seq_len(ncol(z$beta)))mm<-cbind(mm,z$beta[,j]^2)
      colnames(mm)[(ncol(mm)-ncol(z$beta)+1):ncol(mm)]<-paste0("b",seq_len(ncol(z$beta)),"_sq")
      mm<-cbind(mm,s_sq=as.vector(z$s)^2)
      if(prior=="gdp")mm<-cbind(mm,lambda=as.vector(z$hyper))
      if(prior=="lasso")mm<-cbind(mm,L=as.vector(z$hyper))
      mm
    })
    pooled<-do.call(rbind,chains);mu<-colMeans(pooled)
    bm<-do.call(rbind,lapply(chains,function(x)rowsum(x,rep(1:30,each=200))/200))
    mc<-apply(bm,2,sd)/sqrt(nrow(bm))
    target<-unlist(ref$models[[paste(model,prior,sep="_")]]$reference$means)
    for(n in names(target))expect_true(abs(mu[n]-target[n])<6*mc[n]+2e-4,
      info=paste(model,prior,n,"six batch-MCSE plus 0.0002 quadrature allowance"))
  }
})

test_that("GDP complete transition preserves independent reference moments", {
  skip_if_not_installed("jsonlite")
  data<-jsonlite::fromJSON(test_path("fixtures/reference_data.json"))
  ref<-jsonlite::fromJSON(test_path("fixtures/reference_quadrature.json"),simplifyVector=FALSE)
  for(d in 1:2) {
    bs<-as.matrix(read.csv(test_path(paste0("fixtures/gdp_stationary_states_d",d,".csv"))))
    # Independent reference states; no production posterior chains.
    H<-if(d==1)matrix(data$x,ncol=1)else kronecker(diag(2),matrix(data$x,ncol=1))
    y<-if(d==1)data$yr else c(data$yr,data$yi);P<-data$parameters
    set.seed(26000+d)
    z<-cplxshrink:::gdp_invariant_cpp(y,H,d,if(d==2)data$rho_re else 0,
      if(d==2)data$rho_im else 0,bs,P$alpha,P$eta,P$a_s,P$b_s)
    target<-unlist(ref$models[[if(d==1)"real_gdp"else"complex_gdp"]]$reference$means)
    for(j in 1:d)expect_lt(abs(mean(z[,j])-target[paste0("b",j)]),6*sd(z[,j])/sqrt(nrow(z))+2e-4)
    expect_lt(abs(mean(z[,d+1])-target['s']),6*sd(z[,d+1])/sqrt(nrow(z))+2e-4)
  }
})
