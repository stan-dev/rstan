llsurface<- function(dat, ctm, range=2, nsamples=100000){
  
  fit=ctStanFit(datalong = dat,estonly=TRUE,
    ctstanmodel = ctm,optimize=TRUE,#nopriors=TRUE,
    optimcontrol = list(estonly=F,deoptim=F,isloops=0,finishsamples=500),verbose=0)
  
  np=get_num_upars(fit$stanfit$stanfit)
  
  grid=matrix(rnorm(nsamples*np,0,range),ncol=np)
  est=rep(NA,nrow(grid))
  
  for(ri in 1:nrow(grid)){
    est[ri] = log_prob(fit$stanfit$stanfit,grid[ri,])
  }
  
  
  col=1:ncol(grid)
  colt=adjustcolor(col = col,alpha.f = .05)
  matplot(log(-est),grid[,2,drop=F],pch=1,col=colt,xlab= '- log log likelihood', ylab='par value')
  legend('topright',fit$setup$popsetup$parname[match(fit$setup$popsetup$param,1:ncol(grid))],text.col=col,bty='n')
  
  # of=optimizing(object = fit$stanmodel,data= fit$standata,save_iterations=T, hessian=FALSE, iter=1e6, as_vector=FALSE,draws=0,constrained=FALSE,
  #     tol_obj=tol, tol_rel_obj=0,init_alpha=.00001, tol_grad=0,tol_rel_grad=0,tol_param=0,history_size=50)
  # of$par$DRIFT
  # unconstrain_pars(fit$stanfit$stanfit, of$par)
}
