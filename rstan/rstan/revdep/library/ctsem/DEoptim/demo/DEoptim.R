demo.DEoptim <- function(){

  'print.comments' <- function(str){
    star <- "**********"
    cat(paste("\n",star,"\n",str,"\n",star,"\n",sep=""))
  }
  
  'wait' <- function(){
    t <- readline("\nPlease 'q' to quit the demo or any other key to continue...\n")
    if (t == "q") TRUE else FALSE
  }
    
  'Rosenbrock' <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }

  'Wild' <- function(x)
    10 * sin(0.3*x) * sin(1.3*x^2) + 
      0.00001 * x^4 + 0.2 * x + 80

  'demo.1' <- function(){
    r <- DEoptim(Rosenbrock, rep(-10,2), rep(10,2))
    summary(r)
  }

  'demo.2' <- function(){
    r <- DEoptim(Rosenbrock, rep(-10,2), rep(10,2), 
                 control = list(NP = 100, trace = 1))
    summary(r)
  }

  'demo.3' <- function(){
    r <- DEoptim(Rosenbrock, rep(-10,2), rep(10,2), 
                 control = list(NP = 50, itermax = 300, F = 1.5, 
                   CR = 0.2, trace = 1))
    summary(r)
    plot(r, type = 'b')
  }

  'demo.4' <- function(){
    r <- DEoptim(Wild, lower = -50, upper = 50,
                 control = list(NP = 50, trace = 1))
    par(mfrow = c(2,1))
    plot(r, type = 'b')
    plot(r, plot.type = "bestvalit", type = 'l')
  }

  'demo.5' <- function(){
    r <- DEoptim(Wild, lower = -50, upper = 50,
                 control = list(NP = 50, trace = 1, digits = 8))
  }

  str.stop <- "end of the demo"
  tstr <- "\nRun the optimization process for the 'Rosenbrock'"
  tstr <- paste(tstr, "\nBanana function. Search space [-10,10]^2.\n", sep = "")
  print.comments(tstr)
  print(Rosenbrock)
  print(demo.1)
  if (wait()) stop(str.stop) else demo.1()
  
  tstr <- "\nDecrease to 100 the members in the population.\n"
  print.comments(tstr)
  print(demo.2)
  if (wait()) stop(str.stop) else demo.2()
  
  tstr <- "\nIncrease the number of iterations to 300, and"
  tstr <- paste(tstr, "\nmodify crossover and F parameters.\n", sep = "")
  tsts <- paste(tstr, "the result")
  print.comments(tstr)
  print(demo.3)
  if (wait()) stop(str.stop) else demo.3()
  
  tstr <- "\nRun the optimization process for the 'Wild' function."
  tstr <- paste(tstr, "\nSearch space [-50,50].\n", sep = "")
  print.comments(tstr)
  print(Wild)
  plot(Wild, -50, 50, n = 1000,
       main = "DEoptim minimizing 'Wild function'")
  if (wait()) stop(str.stop) else demo.4()

#  tstr <- "\nIncrease the number of printed digits"
#  print.comments(tstr)
#  if (wait()) stop(str.stop) else demo.5()
  
  cat("\n",str.stop,"\n")
}

demo.DEoptim()
  
