test_that("stan_rdump works", {
  skip("Backwards compatibility")

  a <- 1000000000
  b <- rnorm(3)
  c <- rexp(10)
  h <- 1.1e10
  i <- 5.0
  dumpf <- "standumpabc.Rdump"
  stan_rdump(c("a", "b", "c", "h", "i"), file = dumpf)
  text <- paste(readLines(dumpf), collapse = "|")
  expect_true(grepl("a.*<-.*1000000000\\|", text, perl = TRUE))
  expect_true(grepl("h.*<-.*1.1[Ee]\\+*10\\|", text, perl = TRUE))
  on.exit(unlink(dumpf), add = TRUE)
})

test_that("stan_rdump works take 2", {
  skip("Backwards compatibility")

  dat_list <- list(a = array(0, dim = c(3, 2, 1, 0)), d = double(),
                   i = integer(), m = matrix(0, nrow = 3, ncol = 0))
  dumpf <- 'standumpabc.Rdump'
  stan_rdump(names(dat_list), file = dumpf,
                    envir = as.environment(dat_list))
  x <- read_rdump(dumpf)
  on.exit(unlink(dumpf), add = TRUE)

  expect_true(!is.null(x$a))
  expect_true(!is.null(x$d))
  expect_true(!is.null(x$i))
  expect_true(!is.null(x$m))

  expect_true(is.integer(x$a))
  expect_true(is.integer(x$d))
  expect_true(is.integer(x$i))
  expect_true(is.integer(x$m))

  expect_true(length(x$a) == 0)
  expect_true(length(x$d) == 0)
  expect_true(length(x$i) == 0)
  expect_true(length(x$m) == 0)

  expect_equal(dim(x$a), c(3L, 2L, 1L, 0L))
  expect_equal(dim(x$m), c(3L, 0L))
})
