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


test_stan_rdump2 <- function() {
  dat_list <- list(a = array(0, dim = c(3, 2, 1, 0)), d = double(), 
                   i = integer(), m = matrix(0, nrow = 3, ncol = 0))
  dumpf <- 'standumpabc.Rdump'
  rstan::stan_rdump(names(dat_list), file = dumpf, envir = as.environment(dat_list))
  x <- rstan::read_rdump(dumpf)
  unlink(dumpf)
  checkTrue(!is.null(x$a))  
  checkTrue(!is.null(x$d))
  checkTrue(!is.null(x$i))
  checkTrue(!is.null(x$m))

  checkTrue(is.integer(x$a))
  checkTrue(is.integer(x$d))
  checkTrue(is.integer(x$i))
  checkTrue(is.integer(x$m))

  checkTrue(length(x$a) == 0)
  checkTrue(length(x$d) == 0)
  checkTrue(length(x$i) == 0)
  checkTrue(length(x$m) == 0)
  checkEquals(dim(x$a), c(3L, 2L, 1L, 0L))
  checkEquals(dim(x$m), c(3L, 0L))
}

