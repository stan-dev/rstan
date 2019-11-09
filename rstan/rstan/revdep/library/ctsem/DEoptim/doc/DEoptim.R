### R code from vignette source 'DEoptim.Rnw'

###################################################
### code chunk number 1: opt
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: rast
###################################################
rastrigin <- function(x) 
  10*length(x)+sum(x^2-10*cos(2*pi*x))


###################################################
### code chunk number 3: figure1-code
###################################################
library("colorspace")
library("grid")
library("lattice")
jpeg("Rastrigin1.jpg")
x <- y <- seq(-5,5,by=.1)
z <- matrix(nrow=length(x), ncol=length(y))
for(i in 1:length(x)) {
  for(j in 1:length(y)) 
    z[i,j] <- rastrigin(c(x[i],y[j]))
}
xx <- list(fontsize=list(text=15,points=10),
           par.xlab.text=list(cex=2),
           par.ylab.text=list(cex=2),
           axis.text=list(cex=2),
           par.main.text=list(cex=2),
           layout.widths=list(left.padding=.1, right.padding=.1,
             between=0),
           layout.heights=list(top.padding=.1, bottom.padding=.1,
             between=0)
           )
levelplot(z, row.values=x,column.values=y,
          col.regions=sequential_hcl(300), xlab=expression(x[1]), 
          ylab=expression(x[2]),
          par.settings=xx,
          panel=function(z,row.values,column.values,...){
            panel.levelplot(z,row.values,column.values,...);
            panel.points(0,0,pch=21,col="white",cex=2)})

dev.off()
set.seed(123) 


###################################################
### code chunk number 4: prelim
###################################################
library("DEoptim")


###################################################
### code chunk number 5: opt
###################################################
est.ras <- DEoptim(rastrigin,lower=c(-5,-5),upper=c(5,5),
                   control=list(storepopfrom=1, trace=FALSE))


###################################################
### code chunk number 6: figure2-code
###################################################
pushLayout <- function(nr, nc, name="layout") {
  pushViewport(viewport(layout=grid.layout(nr, nc,
                          just="left", widths=unit(rep(2, nc), "null")),
                        name=name))
  for (i in 1:nr) {
    for (j in 1:nc) {
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=j))
      upViewport()
    }
  }
  upViewport()
}
names.vpPath <- names.viewport <- function(x) x$name

with.vpPath <- with.viewport <- function(data, expr, ...) {
  depth <- if (data$name == "ROOT") 0 else downViewport(names(data))
  result <- eval.parent(substitute(expr))
  upViewport(depth)
  invisible(result)
}

getChildren.viewport <- function(x) x$children  

## end functions for making the plots with lattice 
## specify number of cells to fill and number of rows

n <- 6
nr <- 2
nc <- ceiling(n/nr)
xy <- list(fontsize=list(text=12,points=10),
           par.xlab.text=list(cex=1.5),
           par.ylab.text=list(cex=1.5),
           axis.text=list(cex=1.5),
           par.main.text=list(cex=1.5),
           layout.widths=list(left.padding=.1, right.padding=.1,between=0),
           layout.heights=list(top.padding=.1, bottom.padding=.1,between=0))

jpeg("Rastrigin2.jpg")

grid.newpage()
downViewport(pushLayout(nr, nc))
vpt <- current.vpTree(all=FALSE)
plotat <- c(seq(10,50,by=10),1) 
## something strange with Sweave/grid interaction, gen.1 is getting 
## placed in wrong viewport; 'fixed' by permuting plotat above

for(k in 1:n) {
  i <- plotat[k]
  with(getChildren.viewport(vpt)[[k]],
       print(levelplot(z, row.values=x,column.values=y,
                       xlab=expression(x[1]), 
                       ylab=expression(x[2]), colorkey=FALSE,
                       par.settings=xy,between = list(x = .2),
                       col.regions=sequential_hcl(300),
                       main=paste("Generation",i), 
                       panel=function(z,row.values,column.values,...){
                         panel.levelplot(z,row.values,column.values,...);
                         panel.points(est.ras$member$storepop[[i]],
                                      pch=21,fill="black",col=1,cex=.5);
                         panel.points(0,0,pch=21,col="white",cex=1)}),
             newpage = FALSE))
}
dev.off() 


###################################################
### code chunk number 7: opt
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 8: ban
###################################################
genrose.f <- function(x){
n <- length(x)
fval <- 1.0 + sum (100 * (x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
return(fval)
}


###################################################
### code chunk number 9: ban1
###################################################
n <- 10
ans <- DEoptim(fn=genrose.f, lower=rep(-5, n), upper=rep(5, n),
               control=list(NP=100, itermax=4000,trace=FALSE))


###################################################
### code chunk number 10: ban2
###################################################
ans1 <- optim(par=runif(10,-5,5), fn=genrose.f, method="BFGS",
              control=list(maxit=4000))
              


