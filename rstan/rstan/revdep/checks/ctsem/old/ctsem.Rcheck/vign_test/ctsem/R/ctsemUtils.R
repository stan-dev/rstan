suppressOutput <- function(...,verbose=0){
  if(verbose > 0) return(eval(...)) else return(capture.output(eval(...)))
}

naf <-function(x){
  x[is.na(x)] <- FALSE
  return(x)
}

gridplot <- function(m, maxdim=c(3,3),...){
  d=n2mfrow(dim(m)[length(dim(m))])
  d[d>maxdim] <-maxdim[d>maxdim]
  oldpar<-par()
  par(mfrow=d,mar=c(1.1,1.1,1.1,0),mgp=c(.1,.1,0))
  for(i in 1:dim(m)[length(dim(m))]){
    n=colnames(m)[i]
    if(class(m)=='matrix') plot(m[,i],main=ifelse(is.null(n),i,n),col='red',xlab='',ylab='',...)
    if(class(m)=='array') matplot(m[,,i],main=ifelse(is.null(n),i,n),type='l',xlab='',ylab='',...)
  }
  suppressWarnings(do.call(par,oldpar))
}

# helper function to generate an index matrix, or return unique elements of a matrix
indexMatrix<-function(dimension,symmetrical=FALSE,upper=FALSE,lowerTriangular=FALSE, sep=NULL,starttext=NULL,endtext=NULL,
  unique=FALSE,rowoffset=0,coloffset=0,indices=FALSE,diagonal=TRUE,namesvector=NULL){
  if(is.null(namesvector)) namesvector=1:9999
  if(indices==T) sep<-c(",")
  tempmatrix<-matrix(paste0(starttext,namesvector[1:dimension+rowoffset],sep,rep(namesvector[1:dimension+coloffset],each=dimension),endtext),nrow=dimension,ncol=dimension)
  if(upper==TRUE) tempmatrix<-t(tempmatrix)
  if(symmetrical==TRUE) tempmatrix[col(tempmatrix)>row(tempmatrix)] <-t(tempmatrix)[col(tempmatrix)>row(tempmatrix)]
  if(unique==TRUE && symmetrical==TRUE) tempmatrix<-tempmatrix[lower.tri(tempmatrix,diag=diagonal)]
  if(lowerTriangular==TRUE) tempmatrix[col(tempmatrix) > row(tempmatrix)] <- 0
  if(indices==TRUE){
    tempmatrix<-matrix(c(unlist(strsplit(tempmatrix,","))[seq(1,length(tempmatrix)*2,2)],
      unlist(strsplit(tempmatrix,","))[seq(2,length(tempmatrix)*2,2)]),ncol=2)
  }
  return(tempmatrix)
}


relistarrays <- function(flesh, skeleton){
  skelnames <- names(skeleton)
  skelstruc <- lapply(skeleton,dim)
  count=1
  npars <- length(flesh)
  out <- list()
  for(ni in skelnames){
    if(!is.null(skelstruc[[ni]])){
      out[[ni]] <- array(flesh[count:(count+prod(skelstruc[[ni]]))],dim = skelstruc[[ni]])
      count <- count + prod(skelstruc[[ni]])
    } else {
      out[[ni]] <- flesh[count]
      count <- count + 1
    }
  }
  return(out)
}

#convert id's to ascending numeric 
makeNumericIDs <- function(datalong,idName='id',timeName='time'){
  originalid <- unique(datalong[,idName])
  datalong[,idName] <- match(datalong[,idName],originalid)
  
  datalong <- datalong[order(datalong[,idName],datalong[,timeName]),] #sort by id then time
  
  if(any(is.na(as.numeric(datalong[,idName])))) stop('id column may not contain NA\'s or character strings!')
  return(datalong)
}

#' Right multiply a matrix by its transpose.
#'
#' @param x matrix.
#'
#' @return matrix.
#' @export
#'
#' @examples
#' msquare(t(chol(diag(3,4)+1)))
msquare <- function(x) {
  x %*% t(x)
}


crosscov <- function(a,b){
  da <- a-matrix(colMeans(a),nrow=nrow(a),ncol=ncol(a),byrow=TRUE)
  db <- b-matrix(colMeans(b),nrow=nrow(b),ncol=ncol(b),byrow=TRUE)
  t(da) %*% db / (nrow(a)-1)
  
  # cc <- matrix(NA,nrow=nrow(a))
  # for(i in 1:nrow(a)){
  #   cc[i,] <- da[i,] * db[i,]
  
}



#' ctCollapse
#' Easily collapse an array margin using a specified function.
#' @param inarray Input array of more than one dimension.
#' @param collapsemargin Integers denoting which margins to collapse.
#' @param collapsefunc function to use over the collapsing margin.
#' @param ... additional parameters to pass to collapsefunc.
#' @examples
#' testarray <- array(rnorm(900,2,1),dim=c(100,3,3))
#' ctCollapse(testarray,1,mean)
#' @export
ctCollapse<-function(inarray,collapsemargin,collapsefunc,...){
  indims<-dim(inarray)
  out<-array(plyr::aaply(inarray,(1:length(indims))[-collapsemargin],collapsefunc,...,
    .drop=TRUE),dim=indims[-collapsemargin])
  dimnames(out)=dimnames(inarray)[-collapsemargin]
  return(out)
}

rl<-function(x) { #robust logical - wrap checks likely to return NA's in this
  x[is.na(x)] <- FALSE
  return(x)
}




#' Inverse logit
#' 
#' Maps the stan function so the same code works in R.
#'
#' @param x value to calculate the inverse logit for. 
#'
#' @examples
#' inv_logit(-3)
#' @export
inv_logit<-function(x) {
  exp(x)/(1+exp(x))
}



#' ctDensity
#'
#' Wrapper for base R density function that removes outliers and computes 'reasonable' bandwidth and x and y limits.
#' Used for ctsem density plots.
#' 
#' @param x numeric vector on which to compute density.
#' @param bw either 'auto' or a numeric indicating bandwidth.
#' @param plot logical to indicate whether or not to plot the output.
#' @param ... Further args to density.
#' 
#' @examples
#' y <- ctDensity(exp(rnorm(80)))
#' plot(y$density,xlim=y$xlim,ylim=y$ylim)
#' 
#' #### Compare to base defaults:
#' par(mfrow=c(1,2))
#' y=exp(rnorm(10000))
#' ctdens<-ctDensity(y)
#' plot(ctdens$density, ylim=ctdens$ylim,xlim=ctdens$xlim)
#' plot(density(y))
#' @export

ctDensity<-function(x,bw='auto',plot=FALSE,...){
  xlims=stats::quantile(x,probs=c(.05,.95),na.rm=TRUE)
  sd=sd(xlims)
  xlims[1] = xlims[1] - sd
  xlims[2] = xlims[2] + sd
  # x=x[x>xlims[1]*1.2 & x<xlims[2]*1.2]
  # bw=(max(x)-min(x))^1.2 / length(x)^.4 *.4
  if(bw=='auto') bw=min(sd/100,1e-4)
  
  # xlims=stats::quantile(x,probs=c(.01,.99))
  # mid=mean(c(xlims[2],xlims[1]))
  # xlims[1] = xlims[1] - (mid-xlims[1])/8
  # xlims[2] = xlims[2] + (xlims[2]-mid)/8
  
  out1<-stats::density(x,bw=bw,n=5000,from=xlims[1]-sd,to=xlims[2]+sd)
  ylims=c(0,max(out1$y)*1.1)
  
  if(plot) plot(out1$x, out1$y,type='l', xlim=xlims,ylim=ylims,ylab='Density',xlab='Par. Value',...)
  
  return(list(density=out1,xlim=xlims,ylim=ylims))
}



ctDensityList<-function(x,xlimsindex='all',plot=FALSE,smoothness=1,ylab='Density',
  xlab='Par. Value',colvec='auto',ltyvec='auto',probs=c(.05,.95),
  legend=FALSE, legendargs=list(),...){
  
  if(all(xlimsindex=='all')) xlimsindex <- 1:length(x)
  
  for(i in xlimsindex){
    newxlims=stats::quantile(x[[i]],probs=probs,na.rm=TRUE)
    if(i==1) {
      xlims=newxlims
    } 
    else {
      xlims <- range(c(xlims,newxlims))
    }
  }
  sd=sd(xlims)
  xlims[1] = xlims[1] - sd/2
  xlims[2] = xlims[2] + sd/2
  
  bw=abs(max( 
    min( (sd)/length(x[[1]])^.4,sd/50),
    1e-5)) * smoothness
  
  if(all(colvec=='auto')) colvec=1:length(x)
  if(all(ltyvec=='auto')) ltyvec=1:length(x)
  
  # xlims=stats::quantile(x,probs=c(.01,.99))
  # mid=mean(c(xlims[2],xlims[1]))
  # xlims[1] = xlims[1] - (mid-xlims[1])/8
  # xlims[2] = xlims[2] + (xlims[2]-mid)/8
  # browser()
  denslist<-lapply(1:length(x),function(xi) {
    d=stats::density(x[[xi]],bw=bw,n=5000,from=xlims[1]-sd/2,to=xlims[2]+sd/2,na.rm=TRUE)
    # d$y=d$y/ sum(d$y)/range(d$x)[2]*length(d$y)
    return(d)
  })
  ylims=c(0,max(unlist(lapply(denslist,function(li) max(li$y))))*1.1) * ifelse(legend[1]!=FALSE, 1.2,1)
  
  if(plot) {
    plot(denslist[[1]]$x, denslist[[1]]$y,type='l', xlim=xlims,ylim=ylims,ylab=ylab,xlab=xlab,col=colvec[1],lty=ltyvec[1],...)
    if(length(denslist)>1){
      for(ci in 2:length(denslist)){
        points(denslist[[ci]]$x, denslist[[ci]]$y,type='l', col=colvec[ci],lty=ltyvec[ci],...)
      }
    }
    if(all(legend!=FALSE)) {
      if(is.null(legendargs$col)) legendargs$col = colvec
      if(is.null(legendargs$text.col)) legendargs$text.col = colvec
      if(is.null(legendargs$lty)) legendargs$lty = ltyvec
      if(is.null(legendargs$x)) legendargs$x='topright'
      if(is.null(legendargs$bty)) legendargs$bty='n'
      legendargs$legend = legend
      do.call(graphics::legend,legendargs)
    }
  }
  
  return(list(density=denslist,xlim=xlims,ylim=ylims))
}


#' Plots uncertainty bands with shading
#'
#' @param x x values
#' @param y y values
#' @param ylow lower limits of y
#' @param yhigh upper limits of y
#' @param steps number of polygons to overlay - higher integers lead to 
#' smoother changes in transparency between y and yhigh / ylow.
#' @param ... arguments to pass to polygon()
#'
#' @return Nothing. Adds a polygon to existing plot.
#' @export
#'
#' @examples
#' plot(0:100,sqrt(0:100),type='l')
#' ctPoly(x=0:100, y=sqrt(0:100), 
#' yhigh=sqrt(0:100) - runif(101), 
#' ylow=sqrt(0:100) + runif(101),
#' col=adjustcolor('red',alpha.f=.1))
ctPoly <- function(x,y,ylow,yhigh,steps=20,...){
  for(i in 1:steps){
    tylow= y + (ylow-y)*i/steps
    tyhigh= y + (yhigh-y)*i/steps
    xf <- c(x,x[length(x):1])
    yf <- c(tylow,tyhigh[length(tyhigh):1])
    polygon(xf,yf,border=NA,...)
  }
}



#' ctWideNames
#' sets default column names for wide ctsem datasets. Primarily intended for internal ctsem usage.
#' @param n.manifest number of manifest variables per time point in the data.
#' @param Tpoints Maximum number of discrete time points (waves of data, or measurement occasions) 
#' for an individual in the input data structure.
#' @param n.TDpred number of time dependent predictors in the data structure.
#' @param n.TIpred number of time independent predictors in the data structure.
#' @param manifestNames vector of character strings giving column names of manifest indicator variables
#' @param TDpredNames vector of character strings giving column names of time dependent predictor variables
#' @param TIpredNames vector of character strings giving column names of time independent predictor variables
#' @export

ctWideNames<-function(n.manifest,Tpoints,n.TDpred=0,n.TIpred=0,manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){
  
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  manifestnames<-paste0(manifestNames,"_T",rep(0:(Tpoints-1),each=n.manifest))
  if(n.TDpred > 0 && Tpoints > 1) {
    TDprednames<-paste0(TDpredNames,"_T",rep(0:(Tpoints-1),each=n.TDpred))
  } else {
    TDprednames<-NULL
  }
  if (Tpoints > 1) {
    intervalnames<-paste0("dT",1:(Tpoints-1))
  } else {
    intervalnames <- NULL
  }
  if(n.TIpred>0) TIprednames <- paste0(TIpredNames) else TIprednames <- NULL
  return(c(manifestnames,TDprednames,intervalnames,TIprednames))
}

# generates more complex sequences than seq
cseq <- function(from, to, by){
  temp<-c()
  for(i in from){
    temp<-c(temp,seq(i,to,by))
  }
  temp<-sort(temp)
  return(temp)
}

get_stan_params <- function(object) {
  stopifnot(methods::is(object, "stanfit"))
  params <- grep("context__.vals_r", fixed = TRUE, value = TRUE,
    x = strsplit(rstan::get_cppcode(rstan::get_stanmodel(object)), "\n")[[1]])
  params <- sapply(strsplit(params, "\""), FUN = function(x) x[[2]])
  params <- intersect(params, object@model_pars)
  return(params)
}


get_stan_massmat<-function(fit){
  
  spars<-get_stan_params(fit)
  spars2<-c()
  for(pari in spars){
    spars2<-c(spars2,grep(paste0(pari,'['),names(fit@sim$samples[[1]]),fixed=TRUE))
  }
  
  massmat<-list()
  for(chaini in 1:fit@sim$chains){
    temp<-c()
    for(pari in spars2){
      newval<-stats::cov(cbind(fit@sim$samples[[chaini]][[pari]][(fit@sim$warmup - fit@stan_args[[1]]$control$adapt_term_buffer):fit@sim$warmup]))
      names(newval)<-names(fit@sim$samples[[chaini]])[pari]
      temp<-c(temp,newval)
    }
    massmat[[chaini]]<-temp
  }
  return(massmat)
}


