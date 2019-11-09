#' Runs stan, and plots sampling information while sampling.
#'
#' @param object stan model object
#' @param iter Number of iterations
#' @param chains Number of chains
#' @param ... All the other regular arguments to stan()
#' @export
#' @importFrom shiny runApp 
#' @examples
#' \dontrun{
#' #### example 1 
#' scode <- "
#' parameters {
#'   real y[2]; 
#' } 
#' model {
#'   y[1] ~ normal(0, .5);
#'   y[2] ~ double_exponential(0, 2);
#' } 
#' "
#' sm <- stan_model(model_code = scode)
#' fit1 <- stanWplot(object = sm,iter = 100000,chains=2,cores=1)
#' }

stanWplot <- function(object,iter=2000,chains=4,...){
  
  
  tmpdir=tempdir()
  tmpdir=gsub('\\','/',tmpdir,fixed=TRUE)
  
  windows= Sys.info()[1]=='Windows'
  
  stanplot<-function(chains,seed){
    wd<-  paste0("setwd('",tmpdir,"')")

    writeLines(text=paste0(wd,'
    seed<-',seed,';
    chains<-',chains,';
    iter<-',iter,';
    
    notyet<-TRUE
    while(any(notyet==TRUE)){
      Sys.sleep(1);
      samps<-try(data.table::fread(skip="lp__",
      cmd=paste0("grep -v ^# ",seed,"samples_1.csv")),silent=TRUE)
      if(class(samps)[1] != "try-error" && length(samps) > 0) notyet<-FALSE
    }
    varnames<-colnames(samps);
    shiny::runApp(appDir=list(server=function(input, output,session) {
    
    output$chainPlot <- renderPlot({
    parameter<-input$parameter
    refresh <- input$refresh
    begin<-input$begin
    samps<-list()
    for(chaini in 1:chains) {
      samps[[chaini]]<-try(as.matrix(data.table::fread(select = parameter,skip="lp__",
        cmd=paste0("grep -v ^# ",seed,"samples_",chaini,".csv"))),silent=TRUE)
      if(class(samps[[chaini]])[1]=="try-error" || length(samps[[chaini]]) ==0) samps[[chaini]]=samps[[1]][1,,drop=FALSE]
    }
    
    mini<-min(unlist(lapply(1:chains,function(chaini) samps[[chaini]][begin:length(samps[[chaini]]),parameter])),na.rm=T)
    maxi<-max(unlist(lapply(1:chains,function(chaini) samps[[chaini]][begin:length(samps[[chaini]]),parameter])),na.rm=T)
    lengthi<-max(unlist(lapply(1:chains,function(chaini) length(samps[[chaini]][,parameter]))),na.rm=T) #-1:-begin
    
    plot(begin:length(samps[[1]]),
      samps[[1]][begin:length(samps[[1]]),parameter],
      type="l",xlab="",ylab="",main=parameter,
      log=ifelse(parameter %in% c("stepsize__"),"y",""),
      xlim=c(begin,lengthi),
      ylim=c(mini,maxi)
    )
    
    if(chains > 1) for(chaini in 2:chains){
      points(begin:length(samps[[chaini]]),
      samps[[chaini]][begin:length(samps[[chaini]]),parameter], type="l",xlab="",ylab="",main=parameter,col=chaini)
    }
    grid()
    
    })
    },ui=shiny::fluidPage(
    # Application title
    shiny::titlePanel("stan mid-sampling plots..."),
    shiny::sidebarLayout(
    # Sidebar with a slider input for number of observations
    shiny::sidebarPanel(
    shiny::sliderInput("begin", "Start of range:", min = 1,max=iter,value = 1,step=1), 
    shiny::selectInput("parameter", "Choose a parameter:", choices = varnames),
    shiny::actionButton("refresh", "Refresh sample data")
    ),
    
    # Show a plot of the generated distribution
    shiny::mainPanel(
    shiny::plotOutput("chainPlot")
    )
    ))),
    launch.browser=TRUE)
    quit(save="no")'),con=paste0(tmpdir,"/stanplottemp.R"))
    
    if(windows) system(paste0("Rscript --slave --no-restore -e source(\'",tmpdir,"/stanplottemp.R\')"),wait=FALSE) else
      system(paste0(R.home(component = "home"),"/Rscript --slave --no-restore -e source\\(\\\'",tmpdir,"\\/stanplottemp.R\\\'\\)"),wait=FALSE)
    
  }

stanseed<-floor(as.numeric(Sys.time()))

sample_file<-paste0(tmpdir,'/',stanseed,'samples', ifelse(chains==1,'_1',''),'.csv')

stanplot(chains=chains,seed=stanseed)

out=rstan::sampling(object=object,iter=iter,chains=chains,sample_file=sample_file,,...)

for(chaini in 1:chains) system(paste0("rm ",tmpdir,'/',stanseed,"samples_",chaini,".csv"))
system(paste0('rm ',tmpdir,'/stanplottemp.R'))
return(out)
}
