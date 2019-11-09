stopifnot(require("xts"))
stopifnot(require("microbenchmark"))
# Benchmark [.xts using ISO8601 range on large objects
N <- 2e7
s <- 86400*365.25
x <- .xts(1:N, 1.0*seq(s*20, s*40, length.out = N), tzone = "UTC")

# warmup, in case there's any JIT
for (i in 1:2) {
  x["1999/2001",]
}

profile <- FALSE
if (profile) {
  # Use loop if profiling, so microbenchmark calls aren't included
  Rprof(line.profiling = TRUE)
  for(i in 1:10) {
    x[rng,]
  }
  Rprof(NULL)
  print(srp <- summaryRprof())
} else {
  cat("Subset using ISO-8601 range\n")
  microbenchmark(x["1990",], x["1990/",], x["/2009",],
                 x["1990/1994",], x["1990/1999",], x["1990/2009",], times = 5)
}

cat("Subset using integer vector\n")
i001 <- seq(1, N,   1)
i005 <- seq(1, N,   5)
i010 <- seq(1, N,  10)
i050 <- seq(1, N,  50)
i100 <- seq(1, N, 100)
microbenchmark(x[i001,], x[i005,], x[i010,], x[i050,], x[i100,], times = 5)

cat("Subset using logical vector\n")
l001 <- l005 <- l010 <- l050 <- l100 <- logical(N)
l001[i001] <- TRUE
l005[i005] <- TRUE
l010[i010] <- TRUE
l050[i050] <- TRUE
l100[i100] <- TRUE
microbenchmark(x[l001,], x[l005,], x[l010,], x[l050,], x[l100,], times = 5)

cat("Subset using date-time vector\n")
t001 <- index(x)[i001]
t005 <- index(x)[i005]
t010 <- index(x)[i010]
t050 <- index(x)[i050]
t100 <- index(x)[i100]
microbenchmark(x[t001,], x[t005,], x[t010,], x[t050,], x[t100,], times = 5)
