library(ctsem)
library(testthat)

context("ctRasch") #develop some expectations here!



#anomauth
test_that("ctRasch1", {

invlog=function (x) exp(x)/(1 + exp(x))

    
a <- matrix( "a", nrow=1, ncol=1 )
sqrt.q <- matrix( "sqrt.q", nrow=1, ncol=1 )
b <- matrix( "b", nrow=1, ncol=1 )
I <- 3
beta <- matrix( paste0( "beta", 1:I ), ncol=1 )
beta[2] <- 0
lambda <- matrix( rep( 1, I ), nrow=I, ncol=1 )

d <- cbind(rep(1:10,10),1:100,matrix(rbinom(100*I,size=1,prob=invlog(t(t(matrix(rnorm(100*I),nrow=100))))),ncol=I))
colnames(d) <- c('id','time','Y1','Y2','Y3')
  
  
m <- ctModel( n.latent = 1,
              n.manifest = I,
              MANIFESTMEANS = beta,
              LAMBDA = lambda,
              DRIFT = a,
              DIFFUSION = sqrt.q,
              CINT = b,
              type = "stanct" )
m$pars$indvarying <- FALSE
m$pars$indvarying[ m$pars$matrix %in% 'CINT' ] <- TRUE
m$pars$indvarying[ m$pars$matrix %in% 'T0MEANS' ] <- TRUE
m$manifesttype[]=1
row.number <- which( m$pars$matrix %in% "MANIFESTMEANS" )[3]
m$pars$offset[ row.number ] <- 1
m$pars$meanscale[ row.number ] <- 2
set.seed( 1234 )
start <- Sys.time()
r <- ctStanFit( datalong = d,
                ctstanmodel = m,
                iter = 1000,
                chains = 1,
                cores = 1,
                intoverstates = FALSE,
                stationary = FALSE )
print( runtime <- Sys.time() - start )
summary( r )

})
