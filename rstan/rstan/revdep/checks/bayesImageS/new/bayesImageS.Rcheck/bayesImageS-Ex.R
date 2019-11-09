pkgname <- "bayesImageS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bayesImageS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("getBlocks")
### * getBlocks

flush(stderr()); flush(stdout())

### Name: getBlocks
### Title: Get Blocks of a Graph
### Aliases: getBlocks
### Keywords: spatial

### ** Examples

  #Example 1: split a line into 2 blocks
  getBlocks(mask=c(1,1,1,1,0,0,1,1,0), nblock=2)
  
  #Example 2: split a 4*4 2D graph into 4 blocks in order
  #           to use the chequerboard idea for a neighbourhood structure
  #           corresponding to the second-order Markov random field.
  getBlocks(mask=matrix(1, nrow=4, ncol=4), nblock=4)
  
  #Example 3: split a 3*3*3 3D graph into 8 blocks
  #           in order to use the chequerboard idea for a neighbourhood
  #           structure based on the 18 neighbors definition, where the
  #           neighbors of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  mask <- array(1, dim=rep(3,3))
  getBlocks(mask, nblock=8)



cleanEx()
nameEx("getEdges")
### * getEdges

flush(stderr()); flush(stdout())

### Name: getEdges
### Title: Get Edges of a Graph
### Aliases: getEdges
### Keywords: spatial

### ** Examples

  #Example 1: get all edges of a 1D graph. 
  mask <- c(0,0,rep(1,4),0,1,1,0,0)
  getEdges(mask, neiStruc=2)
  
  #Example 2: get all edges of a 2D graph based on neighbourhood structure
  #           corresponding to the first-order Markov random field.
  mask <- matrix(1 ,nrow=2, ncol=3)
  getEdges(mask, neiStruc=c(2,2,0,0))
  
  #Example 3: get all edges of a 2D graph based on neighbourhood structure
  #           corresponding to the second-order Markov random field.
  mask <- matrix(1 ,nrow=3, ncol=3)
  getEdges(mask, neiStruc=c(2,2,2,2))
  
  #Example 4: get all edges of a 3D graph based on 6 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           N,S,E,W, upper and lower adjacencies. To achieve it, there
  #           are several ways, including the two below.
  mask <- array(1, dim=rep(3,3))
  n61 <- matrix(c(2,2,0,0,
                  0,2,0,0,
                  0,0,0,0), nrow=3, byrow=TRUE)
  n62 <- matrix(c(2,0,0,0,
                  0,2,0,0,
                  2,0,0,0), nrow=3, byrow=TRUE)
  e1 <- getEdges(mask, neiStruc=n61)
  e2 <- getEdges(mask, neiStruc=n62)
  e1 <- e1[order(e1[,1], e1[,2]),]
  e2 <- e2[order(e2[,1], e2[,2]),]
  all(e1==e2)
  
  #Example 5: get all edges of a 3D graph based on 18 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  #           To achieve it, there are several ways, including the one below.
  
  n18 <- matrix(c(2,2,2,2,
                  0,2,2,2,
                  0,0,2,2), nrow=3, byrow=TRUE)  
  mask <- array(1, dim=rep(3,3))
  getEdges(mask, neiStruc=n18)



cleanEx()
nameEx("getNeighbors")
### * getNeighbors

flush(stderr()); flush(stdout())

### Name: getNeighbors
### Title: Get Neighbours of All Vertices of a Graph
### Aliases: getNeighbors
### Keywords: spatial

### ** Examples

  #Example 1: get all neighbours of a 1D graph.
  mask <- c(0,0,rep(1,4),0,1,1,0,0,1,1,1)
  getNeighbors(mask, neiStruc=2)
  
  #Example 2: get all neighbours of a 2D graph based on neighbourhood structure
  #           corresponding to the second-order Markov random field.
  mask <- matrix(1, nrow=2, ncol=3)
  getNeighbors(mask, neiStruc=c(2,2,2,2))
  
  #Example 3: get all neighbours of a 3D graph based on 6 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           N,S,E,W, upper and lower adjacencies. To achieve it, there
  #           are several ways, including the two below.
  mask <- array(1, dim=rep(3,3))
  n61 <- matrix(c(2,2,0,0,
                  0,2,0,0,
                  0,0,0,0), nrow=3, byrow=TRUE)
  n62 <- matrix(c(2,0,0,0,
                  0,2,0,0,
                  2,0,0,0), nrow=3, byrow=TRUE)
  n1 <- getNeighbors(mask, neiStruc=n61)
  n2 <- getNeighbors(mask, neiStruc=n62)
  n1 <- apply(n1, 1, sort)
  n2 <- apply(n2, 1, sort)
  all(n1==n2)
  
  #Example 4: get all neighbours of a 3D graph based on 18 neighbours structure
  #           where the neighbours of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  #           To achieve it, there are several ways, including the one below.
  
  n18 <- matrix(c(2,2,2,2,
                  0,2,2,2,
                  0,0,2,2), nrow=3, byrow=TRUE)  
  mask <- array(1, dim=rep(3,3))
  getNeighbors(mask, neiStruc=n18)



cleanEx()
nameEx("gibbsNorm")
### * gibbsNorm

flush(stderr()); flush(stdout())

### Name: gibbsNorm
### Title: Fit a univariate normal (Gaussian) distribution to the observed
###   data.
### Aliases: gibbsNorm

### ** Examples

y <- rnorm(100,mean=5,sd=2)
res.norm <- gibbsNorm(y, priors=list(mu=0, mu.sd=1e6, sigma=1e-3, sigma.nu=1e-3))
summary(res.norm$mu[501:1000])
summary(res.norm$sigma[501:1000])



cleanEx()
nameEx("mcmcPottsNoData")
### * mcmcPottsNoData

flush(stderr()); flush(stdout())

### Name: mcmcPottsNoData
### Title: Simulate pixel labels using chequerboard Gibbs sampling.
### Aliases: mcmcPottsNoData

### ** Examples

# Swendsen-Wang for a 2x2 lattice
neigh <- matrix(c(5,2,5,3,  1,5,5,4,  5,4,1,5,  3,5,2,5), nrow=4, ncol=4, byrow=TRUE)
blocks <- list(c(1,4), c(2,3))
res.Gibbs <- mcmcPottsNoData(0.7, 3, neigh, blocks, niter=200)
res.Gibbs$z
res.Gibbs$sum[200]



cleanEx()
nameEx("swNoData")
### * swNoData

flush(stderr()); flush(stdout())

### Name: swNoData
### Title: Simulate pixel labels using the Swendsen-Wang algorithm.
### Aliases: swNoData

### ** Examples

# Swendsen-Wang for a 2x2 lattice
neigh <- matrix(c(5,2,5,3,  1,5,5,4,  5,4,1,5,  3,5,2,5), nrow=4, ncol=4, byrow=TRUE)
blocks <- list(c(1,4), c(2,3))
res.sw <- swNoData(0.7, 3, neigh, blocks, niter=200)
res.sw$z
res.sw$sum[200]



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
