
model_code <- "model { \n y ~ normal(0, 1); \n}"
cat(model_code, file = "tmp.stan")

a <- c(1, 3, 5)
b <- matrix(1:10, ncol = 2)
c <- array(1:18, dim = c(2, 3, 3))
dump(c("a", "b", "c"), file = "dumpabc.Rdump")
stan_rdump(c("a", "b", "c"), file = "standumpabc.Rdump")
d <- factor(c("a", "b", "b", "a", "a"))
e <- factor(c(1, 3, 3, 1, 1))
f <- c(TRUE, TRUE, FALSE)
g <- 1000000000
h <- 1.1e+10
stan_rdump(c("d", "e", "f"), file = "standumpabc.Rdump", append = TRUE)
stan_rdump(c("g", "h"), file = "standumpabc.Rdump", append = TRUE)

cc <- c(
  "# comment line 1",
  " no comments line 1",
  "# comment line 2",
  "# comment line 3",
  "# comment line 4",
  "# comment line 5",
  " not comments line 2",
  "# comment line 6",
  "not comments #comment line 7",
  "not comments at the end of file"
)
cat(file = "cc.csv", paste(cc, collapse = "\n"), "\n")
cc2 <- c("# line 1", "#line 2", "# line 3", "#line4", "c,a,b", "1,2,3")
cat(file = "cc2.csv", paste(cc2, collapse = "\n"), "\n")

test_that("get_model_strcode works", {
  model_code <- "model { \n y ~ normal(0, 1); \n}"
  code <- "parameters { real y; } model{y ~ normal(0,1);}"
  # FIXME: Warning - "incomplete final line found on 'tmp.stan'"
  str1 <- suppressWarnings(rstan:::get_model_strcode("tmp.stan"))
  str2 <- rstan:::get_model_strcode(model_code = code)
  str3 <- rstan:::get_model_strcode(model_code = "code")
  str4 <- rstan:::get_model_strcode(model_code = "parameters {real y;} model {y ~ normal(0,1); }")

  mname1 <- attr(str1, "model_name2")
  mname2 <- attr(str2, "model_name2")
  mname3 <- attr(str3, "model_name2")
  mname4 <- attr(str4, "model_name2")
  expect_equal(mname1, "tmp")
  expect_equal(mname2, "code")
  expect_equal(mname3, "code")
  expect_equal(mname4, "anon_model")

  attributes(str1) <- NULL
  attributes(str2) <- NULL
  attributes(str3) <- NULL
  expect_equal(str1, model_code)
  expect_equal(str2, code)
  expect_equal(str3, code)

  model_code <- "model { \n y ~ normal(0, 1); \n}"
  # cat(model_code, file = 'tmp.stan')
  # FIXME: Warning - "incomplete final line found on 'tmp.stan'"
  expect_equal(
    model_code,
    suppressWarnings(rstan:::read_model_from_con("tmp.stan"))
  )
  attr(model_code, "model_name2") <- "tmp"
  # FIXME: Warning - "incomplete final line found on 'tmp.stan'"
  expect_equal(
    model_code,
    suppressWarnings(rstan:::get_model_strcode("tmp.stan"))
  )
  attr(model_code, "model_name2") <- "model_code"
  expect_equal(model_code, rstan:::get_model_strcode(model_code = model_code))
  expect_error(rstan:::get_model_strcode())
})

test_that("is_legal_stan_name works", {
  expect_false(rstan:::is_legal_stan_vname("7dd"))
  expect_false(rstan:::is_legal_stan_vname("model"))
  expect_false(rstan:::is_legal_stan_vname("private"))
  expect_false(rstan:::is_legal_stan_vname("hello__"))
  expect_true(rstan:::is_legal_stan_vname("y"))
})

test_that("is_named_list works", {
  expect_false(rstan:::is_named_list(c(2, 3)))
  expect_false(rstan:::is_named_list(list(3, 4)))
  expect_false(rstan:::is_named_list(list(a = 3, 4)))
  expect_true(rstan:::is_named_list(list(a = 3, b = 4)))
})

test_that("data_preprocess works", {
  lst <- list(
    z = c(1L, 2L, 4L),
    a = 1:100,
    b = matrix(1:9 / 9, ncol = 3),
    c = structure(1:100, .Dim = c(5, 20)),
    g = array(c(3, 3, 9, 3, 3, 4, 5, 6, 9, 8, 0, 2), dim = c(2, 2, 3)),
    d = 1:100 + .1
  )
  lst <- rstan:::data_preprocess(lst)
  lst2 <- lst
  lst2$f <- matrix(c(3, NA, NA, NA, 3, 4), ncol = 3)
  lst3 <- lst
  lst3$h <- gl(3, 4)
  # Warns about `h` being not numeric.
  # FIXME: Should the warning be tested explicitly?
  lst4 <- suppressWarnings(rstan:::data_preprocess(lst3))

  # Keep the dimension information.
  expect_equal(dim(lst$g), c(2, 2, 3))
  # Do as.integer when appropriate.
  expect_true(is.integer(lst$z))
  # Don't do as.integer when it is not appropriate.
  expect_true(is.double(lst$b))
  # Stop if data have NA.
  expect_error(rstan:::data_preprocess(lst2))
  # Check if h is removed.
  expect_named(lst4, c("z", "a", "b", "c", "g", "d"))
})

test_that("data_preprocess works take 2", {
  # A list of array as an element of the data list.
  I <- 3
  J <- 4
  K <- 5
  a <- lapply(1:I, function(i) rnorm(J))
  b <- lapply(1:I, function(i) matrix(rnorm(J * K), ncol = K))
  d <- lapply(1:I, function(i) rnorm(1, i))
  e <- lapply(1:I, function(i) rpois(J, 1) + 1.0)
  lst2 <- rstan:::data_preprocess(list(a = a, b = b, d = d, e = e))
  expect_equal(dim(lst2$a), c(I, J))
  expect_equal(dim(lst2$b), c(I, J, K))
  expect_equal(dim(lst2$d), c(I, 1))
  expect_true(is.integer(lst2$e[1, 1]))
  expect_false(is.integer(e[[1]][1]))

  a1 <- lapply(1:I, function(i) rnorm(J))
  a2 <- lapply(1:I, function(i) rnorm(J))
  expect_error(rstan:::data_preprocess(list(a = list(a1 = a1, a2 = a2))))
})

test_that("read_rdump works", {
  l <- rstan:::read_rdump("dumpabc.Rdump")
  expect_equal(l$a, c(1, 3, 5))
  expect_equal(l$b, matrix(1:10, ncol = 2))
  expect_equal(l$c, array(1:18, dim = c(2, 3, 3)))
})

test_that("stan_rdump works", {
  l <- rstan:::read_rdump("standumpabc.Rdump")
  expect_equal(l$a, c(1, 3, 5))
  expect_equal(l$b, matrix(1:10, ncol = 2))
  expect_equal(l$c, array(1:18, dim = c(2, 3, 3)))
  expect_equal(l$d, c(1, 2, 2, 1, 1))
  expect_equal(l$e, c(1, 2, 2, 1, 1))
  expect_equal(l$f, c(1, 1, 0))
  text <- paste(readLines("standumpabc.Rdump"), collapse = "|")
  expect_true(grepl("1000000000", text))
  expect_equal(l$h, 1.1e10)
})

test_that("seq_array_ind works", {
  a <- rstan:::seq_array_ind(numeric(0))
  expect_length(a, 0)
  # By default, col_major is FALSE.
  b <- rstan:::seq_array_ind(2:5, col_major = TRUE)
  c <- arrayInd(1:prod(2:5), .dim = 2:5)
  expect_equal(b, c)
  d <- rstan:::seq_array_ind(2:3, col_major = FALSE)
  e <- matrix(c(1, 1, 1, 2, 1, 3, 2, 1, 2, 2, 2, 3),
    nrow = 6, byrow = TRUE
  )
  expect_equal(d, as.array(e))
  f <- rstan:::seq_array_ind(1)
  expect_equal(f, array(1, dim = c(1, 1)))
})

test_that("flatnames works", {
  names <- c("alpha", "beta", "gamma", "delta")
  dims <- list(alpha = integer(0), beta = c(2, 3), gamma = c(2, 3, 4),
               delta = c(5))
  fnames <- rstan:::flatnames(names, dims)
  expect_equal(
    fnames,
    c(
      "alpha", "beta[1,1]", "beta[1,2]", "beta[1,3]",
      "beta[2,1]", "beta[2,2]", "beta[2,3]",
      "gamma[1,1,1]", "gamma[1,1,2]", "gamma[1,1,3]", "gamma[1,1,4]",
      "gamma[1,2,1]", "gamma[1,2,2]", "gamma[1,2,3]", "gamma[1,2,4]",
      "gamma[1,3,1]", "gamma[1,3,2]", "gamma[1,3,3]", "gamma[1,3,4]",
      "gamma[2,1,1]", "gamma[2,1,2]", "gamma[2,1,3]", "gamma[2,1,4]",
      "gamma[2,2,1]", "gamma[2,2,2]", "gamma[2,2,3]", "gamma[2,2,4]",
      "gamma[2,3,1]", "gamma[2,3,2]", "gamma[2,3,3]", "gamma[2,3,4]",
      "delta[1]", "delta[2]", "delta[3]", "delta[4]", "delta[5]"
    )
  )
  names2 <- c("alpha")
  dims2 <- list(alpha = integer(0))
  fnames2 <- rstan:::flatnames(names2, dims2)
  expect_equal(fnames2, "alpha")
})

test_that("idx_col2rowm works", {
  d <- integer(0)
  idx <- rstan:::idx_col2rowm(d)
  expect_equal(idx, 1)
  d2 <- 8
  idx2 <- rstan:::idx_col2rowm(d2)
  expect_equal(idx2, 1:8)
  d3 <- c(3, 4, 5)
  idx3 <- rstan:::idx_col2rowm(d3)
  yidx3 <- c(
    1, 13, 25, 37, 49, 4, 16, 28, 40, 52, 7, 19, 31, 43, 55,
    10, 22, 34, 46, 58, 2, 14, 26, 38, 50, 5, 17, 29, 41, 53, 8, 20, 32, 44, 56,
    11, 23, 35, 47, 59, 3, 15, 27, 39, 51, 6, 18, 30, 42, 54, 9, 21, 33, 45, 57,
    12, 24, 36, 48, 60
  )
  expect_equal(idx3, yidx3)
})

test_that("idx_row2colm works", {
  d <- integer(0)
  idx <- rstan:::idx_row2colm(d)
  expect_equal(idx, 1)
  d2 <- 8
  idx2 <- rstan:::idx_row2colm(d2)
  expect_equal(idx2, 1:8)
  d3 <- c(3, 4, 5)
  idx3 <- rstan:::idx_row2colm(d3)
  yidx3 <- c(
    1, 21, 41, 6, 26, 46, 11, 31, 51, 16, 36, 56, 2, 22, 42, 7, 27, 47, 12, 32, 52,
    17, 37, 57, 3, 23, 43, 8, 28, 48, 13, 33, 53, 18, 38, 58, 4, 24, 44, 9, 29, 49,
    14, 34, 54, 19, 39, 59, 5, 25, 45, 10, 30, 50, 15, 35, 55, 20, 40, 60
  )
  expect_equal(idx3, yidx3)
})

test_that("pars_total_indexes works", {
  names <- "alpha0"
  dims <- list(alpha0 = c(2, 3))
  fnames <- rstan:::flatnames(names, dims)
  tidx <- rstan:::pars_total_indexes(names, dims, fnames, "alpha0")
  tidx_attr1 <- attr(tidx[[1]], "row_major_idx")
  attributes(tidx[[1]]) <- NULL
  expect_equal(unname(tidx[[1]]), 1:6)
  expect_equal(unname(tidx_attr1), c(1, 3, 5, 2, 4, 6))
  names2 <- c(names, "alpha")
  dims2 <- c(dims, list(alpha = 8))
  fnames2 <- rstan:::flatnames(names2, dims2)
  tidx2 <- rstan:::pars_total_indexes(names2, dims2, fnames2, "alpha")
  tidx2_attr1 <- attr(tidx2[[1]], "row_major_idx")
  attributes(tidx2[[1]]) <- NULL
  expect_equal(unname(tidx2[[1]]), 6 + 1:8)
  expect_equal(unname(tidx2_attr1), 6 + 1:8)
  names3 <- c(names2, "p")
  dims3 <- c(dims2, list(p = integer(0)))
  fnames3 <- rstan:::flatnames(names3, dims3)
  tidx3 <- rstan:::pars_total_indexes(names3, dims3, fnames3, "p")
  tidx3_attr1 <- attr(tidx3[[1]], "row_major_idx")
  attributes(tidx3[[1]]) <- NULL
  expect_equal(unname(tidx3[[1]]), 15)
  expect_equal(unname(tidx3_attr1), 15)
})

test_that("multi_idx_row2colm works", {
  expect_equal(rstan:::multi_idx_row2colm(list(integer(0))), 1)
  dims <- list(c(3), c(2, 3), integer(0), c(2))
  col_idx <- rstan:::multi_idx_row2colm(dims)
  target <- c(1, 2, 3, 4, 7, 5, 8, 6, 9, 10, 11, 12)
  expect_equal(col_idx, target)

  fnames <- c(
    "alpha[1]", "alpha[2]", "alpha[3]",
    "alpha2[1,1]", "alpha2[1,2]", "alpha2[1,3]",
    "alpha2[2,1]", "alpha2[2,2]", "alpha2[2,3]",
    "p", "theta[1]", "theta[2]"
  )
  fnames_colm <- c(
    "alpha[1]", "alpha[2]", "alpha[3]",
    "alpha2[1,1]", "alpha2[2,1]", "alpha2[1,2]", "alpha2[2,2]",
    "alpha2[1,3]", "alpha2[2,3]", "p", "theta[1]", "theta[2]"
  )
  expect_equal(fnames[col_idx], fnames_colm)
})

test_that("mklist works", {
  x <- 3:5
  y <- array(1:9, dim = c(3, 3))
  z <- list(p = 3)
  f <- function() {
    TRUE
  }
  a <- list(x = x, y = y)
  b <- rstan:::mklist(c("x", "y"))
  expect_identical(a, b)
  c <- list(x = x, y = y, z = z)
  d <- rstan:::mklist(c("x", "y", "z"))
  expect_identical(c, d)
  expect_error(rstan:::mklist(c("x", "f")))

  p <- 4
  fun <- function() {
    p <- 3
    rstan:::mklist("p")
  }
  expect_equal(fun()$p, 3, ignore_attr = "names")
  expect_equal(rstan:::mklist("p")$p, 4, ignore_attr = "names")
})

# test_makeconf_path <- function() {
#  p <- makeconf_path()
#  checkTrue(file.exists(makeconf_path()))
# }

test_that("config_argss works", {
  # (chains, iter, warmup, thin, init, seed, sample_file, ...)
  a <- rstan:::config_argss(3, 100, 10, 3, 0, 0, "a.csv",
    algorithm = "NUTS",
    control = NULL,
    chain_id = 4
  )
  expect_length(a, 3)
  expect_equal(a[[1]]$init, "0")
  expect_equal(a[[1]]$chain_id, 4)
  expect_equal(a[[3]]$chain_id, 6)
  b <- rstan:::config_argss(3, 100, 10, 3, "0", 10, "a.csv",
    algorithm = "NUTS",
    control = NULL
  )
  expect_equal(b[[3]]$chain_id, 3)
  expect_equal(b[[1]]$init, "0")
  c <- rstan:::config_argss(3, 100, 10, 3, "random", 10, "a.csv",
    algorithm = "HMC",
    control = list(adapt_engaged = FALSE)
  )
  expect_equal(c[[1]]$init, "random")
  d <- rstan:::config_argss(4, 100, 10, 3, "random", 10, "a.csv",
    chain_id = c(3, 2, 1), algorithm = "Metropolis",
    control = list(adapt_engaged = FALSE)
  )
  expect_equal(d[[3]]$chain_id, 1)
  expect_equal(d[[4]]$chain_id, 4)
  expect_error(rstan:::config_argss(3, 100, 10, 3, "random", 10, NA,
    algorithm = "NUTS", control = NULL, chain_id = c(3, 3)
  ))
  b <- rstan:::config_argss(3, 100, 10, 3, 0, "12345", "a.csv",
    algorithm = "NUTS", chain_id = 4, control = NULL
  )
  expect_equal(b[[1]]$seed, "12345")
  expect_error(rstan:::config_argss(3, 100, 10, 3, 0, "a12345", "a.csv",
    algorithm = "NUTS", control = NULL, chain_id = 4
  ))
  expect_error(rstan:::config_argss(3, 100, 10, 3, 0, "1a2345", "a.csv",
    algorithm = "NUTS", control = NULL,
    chain_id = 4
  ))
})

test_that("data_list2array works", {
  d <- list(y = rnorm(20))
  d2 <- rstan:::data_list2array(d)
  expect_equal(d2, array(d$y, dim = c(1, 20)))

  I <- 4
  J <- 5
  K <- 6
  b <- lapply(1:I, function(i) rnorm(J))
  b2 <- rstan:::data_list2array(b)
  b3 <- data.matrix(do.call(rbind, b))
  expect_equal(b2, b3)

  a <- lapply(1:I, function(i) array(rnorm(J * K), dim = c(J, K)))
  a2 <- rstan:::data_list2array(a)
  for (i in 1:I) {
    expect_equal(a[[i]], a2[i, , ])
  }
  expect_equal(a[[4]][5, 6], a2[4, 5, 6])
  expect_equal(a[[3]][4, 2], a2[3, 4, 2])
  expect_equal(a[[2]][1, 6], a2[2, 1, 6])
})

test_that("read_comments works", {
  a1 <- rstan:::read_comments("cc.csv", 5L)
  expect_equal(a1[1], "# comment line 1")
  expect_equal(a1[5], "# comment line 5")
  a2 <- rstan:::read_comments("cc.csv", 3L)
  expect_length(a2, 3L)
  a3 <- rstan:::read_comments("cc.csv", -1)
  expect_length(a3, 7L)
})

test_that("read_csv_header", {
  h1 <- rstan:::read_csv_header("cc2.csv")
  expect_equal(attr(h1, "lineno"), 5)
  expect_equal(h1, "c,a,b", ignore_attr = TRUE)
})

test_that("get_dims_from_fnames works", {
  names <- c("alpha", "beta2", "g2amma", "theta0")
  dims <- list(c(2L, 3L), integer(0L), c(4L), c(3L, 5L, 4L))
  fnames <- rstan:::flatnames(names, dims)
  fnames_d <- rstan:::sqrfnames_to_dotfnames(fnames)
  unames <- rstan:::unique_par(fnames_d)
  dims2 <- lapply(
    unames,
    function(n) {
      fnames_d2 <- fnames_d[sapply(fnames_d, function(i) grepl(n, i))]
      # The above line works here since all parameters are not nested.
      # it would be problematic if say we have another parameter `p`,
      # since p is also part of `alpha`.
      rstan:::get_dims_from_fnames(fnames_d2)
    }
  )
  expect_equal(dims, dims2)
})

test_that("Converting from/to dotfnames/sqrfnames works", {
  dn <- c("alpha", "beta.1", "beta.2", "gamma.1.2", "gamma.1.4")
  sn <- rstan:::dotfnames_to_sqrfnames(dn)
  dn2 <- rstan:::sqrfnames_to_dotfnames(sn)
  expect_equal(dn, dn2)
})

test_that("Can convert parameter from vector to list", {
  v <- c(2.3, 3.4, 4.5, (1:8) / 9, 3.1415)
  pars <- c("alpha", "beta", "gamma", "delta")
  dims <- list(integer(0), c(2), c(2, 4), 1)
  vl <- rstan:::rstan_relist(v, rstan:::create_skeleton(pars, dims))
  alpha <- 2.3
  beta <- array(v[2:3], dim = 2)
  gamma <- array(v[4:11], dim = c(2, 4))
  delta <- array(v[12], dim = 1)
  expect_length(vl, 4)
  expect_equal(vl[[1]], alpha)
  expect_equal(vl[[2]], beta)
  expect_equal(vl[[3]], gamma)
  expect_equal(vl[[4]], delta)
})

test_that("remove_empty_pars works", {
  pars <- c("alpha", "beta", "gamma", "eta", "xi")
  dims <- list(integer(0), c(2), c(2, 4), 0, c(2, 0))
  names(dims) <- pars
  expect_equal(rstan:::remove_empty_pars(pars[1:2], dims), pars[1:2])
  expect_equal(rstan:::remove_empty_pars(pars[1:4], dims), pars[1:3])
  expect_equal(rstan:::remove_empty_pars("beta[1]", dims), "beta[1]")
  expect_equal(rstan:::remove_empty_pars("eta", dims), character(0))
})

unlink("tmp.stan")
unlink("dumpabc.Rdump")
unlink("standumpabc.Rdump")
unlink("cc.csv")
unlink("cc2.csv")
