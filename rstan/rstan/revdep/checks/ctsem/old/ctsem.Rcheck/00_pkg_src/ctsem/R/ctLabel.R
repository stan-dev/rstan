#Function used internally by ctsem to create paramater names for matrices

ctLabel<-function(matrixname, n.latent, n.manifest, n.TDpred, n.TIpred, Tpoints, manifestNames, latentNames, TDpredNames, TIpredNames){
  
  if(matrixname=="T0MEANS") out <- matrix(paste0("T0mean_",latentNames[1:n.latent]),ncol=1)
  if(matrixname=="T0VAR") out <- indexMatrix(dimension=n.latent,starttext="T0var_",lowerTriangular=TRUE,sep="_",namesvector=latentNames)
  if(matrixname=='LAMBDA') out <- matrix(paste0('lambda_',rep(1:n.manifest,times=n.latent),rep(latentNames[1:n.latent],each=n.manifest)),nrow=n.manifest)
  if(matrixname=="MANIFESTMEANS") out <- matrix(paste0('manifestmeans_',manifestNames[1:n.manifest]),nrow=n.manifest,ncol=1)
  if(matrixname=="MANIFESTVAR") out <- indexMatrix(dimension=n.manifest,starttext="manifestvar_",lowerTriangular=TRUE,sep="_",namesvector=manifestNames)
  if(matrixname=="DRIFT") out <- indexMatrix(dimension=n.latent,starttext="drift_",symmetrical=FALSE,sep="_",namesvector=latentNames)
  if(matrixname=="CINT") out <- matrix(paste0("cint_",latentNames[1:n.latent]),ncol=1)
  if(matrixname=="DIFFUSION") out <- indexMatrix(dimension=n.latent,starttext="diffusion_",lowerTriangular=TRUE,sep="_",namesvector=latentNames)
  if(matrixname=="TRAITVAR") out <- indexMatrix(dimension=n.latent,starttext="traitvar_",lowerTriangular=TRUE,sep="_",namesvector=latentNames)
  if(matrixname=="T0TRAITEFFECT") out <- indexMatrix(dimension=n.latent,starttext="T0traiteffect_",symmetrical=FALSE,sep="_",namesvector=latentNames,endtext='Trait')
  if(matrixname=="MANIFESTTRAITVAR") out <- indexMatrix(dimension=n.manifest,starttext="manifesttraitvar_",lowerTriangular=TRUE,sep="_",namesvector=manifestNames)
  if(matrixname=="MANIFESTMEANS") out <- matrix(paste0("manifestmeans_",manifestNames[1:n.manifest]),ncol=1)
  
  if(matrixname=="TDPREDEFFECT") out <- matrix(paste0("TDeffect_",latentNames[1:n.latent],"_",rep(TDpredNames[1:n.TDpred],each=n.latent)),nrow=n.latent,ncol=n.TDpred)
  
  if(matrixname=="TDPREDMEANS") out <- matrix(paste0("mean_",rep(TDpredNames[1:n.TDpred],each=Tpoints),"T",0:(Tpoints-1)),ncol=1)
   
  if(matrixname=="TIPREDMEANS") out <- matrix(paste0("mean_",TIpredNames),ncol=1)

  if(matrixname=="T0TDPREDCOV") out <- matrix(paste0(
    "T0TDPREDCOV_",
    latentNames[1:n.latent],
    "_",
    rep(TDpredNames[1:n.TDpred],each=n.latent*(Tpoints)),"_T",
    rep(0:(Tpoints-1),each=n.latent)),nrow=n.latent,ncol=n.TDpred*(Tpoints))
  
  if(matrixname=="TIPREDEFFECT") out <- matrix(paste0("TIeffect_",
    latentNames[1:n.latent],"_",rep(TIpredNames[1:n.TIpred],each=n.latent)),ncol=n.TIpred,nrow=n.latent)
  
  if(matrixname=="T0TIPREDEFFECT") out <- matrix(paste0(  
    "T0TIeffect_",
    latentNames[1:n.latent],
    "_",
    rep(TIpredNames[1:n.TIpred],each=n.latent)),nrow=n.latent,ncol=n.TIpred)
  
  
  if(matrixname=="TDPREDVAR")  {
    
    out <- matrix(paste0(
    rep(TDpredNames,each=Tpoints-1), #row predictor
    "T", 
    0:(Tpoints-1), #row time
    "_", 
    rep(TDpredNames,each=((Tpoints)*n.TDpred*(Tpoints-1))), #col predictor
    "T", 
    rep(0:(Tpoints-1),each=(Tpoints)*n.TDpred), #col time
    "_cov"),
    nrow=(n.TDpred*(Tpoints)), ncol=(n.TDpred*(Tpoints)))
  
    out[upper.tri(out)]<- 0 #ensure lower triangular
    
  }
  
  
  
  if(matrixname=="TIPREDVAR") {
    out <- matrix(paste0(
    TIpredNames, 
    "_", 
    rep(TIpredNames,each=n.TIpred), 
    "_cov"),
    nrow=n.TIpred, ncol=n.TIpred)
    
    out[upper.tri(out)]<- 0 #ensure lower triangular
  }
  
  
  if(matrixname=="TDTIPREDCOV") out <- matrix(paste0(
    rep(TDpredNames,each=Tpoints-1), #row predictor
    "T",
    0:(Tpoints-1), #row time
    "_",
    rep(TIpredNames,each=n.TDpred*(Tpoints)),
    "_cov"),
    nrow=(n.TDpred*(Tpoints)), ncol=n.TIpred)
  
  
  



if(matrixname=="TRAITTDPREDCOV"){
  if(n.TDpred>0) {
    TRAITTDPREDCOV<-matrix(paste0("traitTDcov_",latentNames[1:n.latent],"_",
      rep(TDpredNames[1:n.TDpred],each=(Tpoints)*n.latent),
      paste0("T",rep(0:(Tpoints-1),each=n.latent))),
      nrow=n.latent,ncol=n.TDpred*(Tpoints))
    out <-TRAITTDPREDCOV
  }
}  
  

  
  
  

return(out)
}
