test_stan_rdump <- function() {
  a <- 1000000000
  b <- rnorm(3)
  c <- rexp(10)
  h <- 1.1e10
  i <- 5.0
  dumpf <- 'standumpabc.Rdump'
  rstan:::stan_rdump(c("a", "b", "c", "h", "i"), file = dumpf) 
  text <- paste(readLines(dumpf), collapse = "|")
  checkTrue(grepl("a.*<-.*1000000000\\|", text, perl = TRUE))
  checkTrue(grepl("h.*<-.*1.1[Ee]\\+*10\\|", text, perl = TRUE))
  unlink(dumpf)
} 

