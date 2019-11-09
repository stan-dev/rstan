# This benchmark script is based on one written by Dirk Eddelbuettel 
# while he was porting DEoptim to Rcpp.  His efforts in code review
# and benchmarking are greatly appreciated.
demo.DEbenchmark <- function() {
    Rosenbrock <- function(x){
        x1 <- x[1]
        x2 <- x[2]
        100 * (x2 - x1 * x1)^2 + (1 - x1)^2
    }
    
    Wild <- function(x) {       ## 'Wild' function, global minimum at about -15.81515
        sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80)/length(x)
    }
    
    Rastrigin <- function(x) {
        sum(x+2 - 10 * cos(2*pi*x)) + 20
    }
    
    Genrose <- function(x) {        ## One generalization of the Rosenbrock banana valley function (n parameters)
        n <- length(x)
        1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
    }

    maxIt <- 250                        # not excessive but so that we get some run-time on simple problems
    
    basicDE <- function(n, maxIt, fun) DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                control=list(NP=10*n, itermax=maxIt, trace=FALSE))#, bs=TRUE))
    adaptDE <- function(n, maxIt, fun) DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                control=list(NP=10*n, itermax=maxIt, trace=FALSE, strategy=6))#, bs=TRUE))
    
    runPair <- function(n, maxIt, fun) {
        
        gc()
        set.seed(42)
        bt <- mean(replicate(10, system.time(invisible(basicDE(n, maxIt, fun)))[3]), trim=0.1)
        
        gc()
        set.seed(42)
        ct <- mean(replicate(10, system.time(invisible(adaptDE(n, maxIt, fun)))[3]), trim=0.1)
        
        return(data.frame(DE=bt, JADE=ct))
    }
    
    cat("# At", format(Sys.time()), "\n")
    
    reps <- c(5, 10, 20)
    
    res <- rbind(
            do.call(rbind, lapply(reps, runPair, maxIt, function(...) Rosenbrock(...))),
            do.call(rbind, lapply(reps, runPair, maxIt, function(...) Rastrigin(...))),
            do.call(rbind, lapply(reps, runPair, maxIt, function(...) Wild(...))),
            do.call(rbind, lapply(reps, runPair, maxIt, function(...) Genrose(...)))
    )
    res <- rbind(res, colMeans(res))
    rownames(res) <- c( paste("Rosenbrock", reps, sep=""),
                        paste("Rastrigin", reps, sep=""),
                        paste("Wild", reps, sep=""),
                        paste("Genrose", reps, sep=""),
                        "MEANS")
    res
}

demo.DEbenchmark()