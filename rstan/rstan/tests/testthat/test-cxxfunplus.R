test_that("cxxfunctionplus propagates compilation errors", {
  # https://github.com/stan-dev/rstan/issues/1197
  local_mocked_bindings(
    cxxfunction = function(...) stop("compilation failed"),
    rstan_options = function(...) FALSE
  )
  n_sinks <- sink.number()
  expect_error(cxxfunctionplus(), "compilation failed")
  expect_equal(sink.number(), n_sinks)
})

test_that("cxxfunctionplus cleans up after a successful compile", {
  # https://github.com/stan-dev/rstan/issues/1199
  skip_if(Sys.getenv("USE_CXX17") != "")
  local_mocked_bindings(
    cxxfunction = function(...) "fx",
    dso_path = function(...) "fake.so",
    get_CXX = function(...) "g++",
    get_makefile_flags = function(...) "",
    rstan_options = function(...) FALSE
  )
  n_sinks <- sink.number()
  expect_s4_class(cxxfunctionplus(module_name = ""), "cxxdso")
  expect_equal(Sys.getenv("USE_CXX17"), "")
  expect_equal(sink.number(), n_sinks)
})

test_that("cxxfunctionplus skips user Makevars with -march=native on Windows", {
  # https://github.com/stan-dev/rstan/issues/1199
  skip_on_os(c("mac", "linux", "solaris"))
  makevars_during <- NULL
  makevars_size <- NULL
  local_mocked_bindings(
    .warn_march_makevars = function() TRUE,
    cxxfunction = function(...) {
      makevars_during <<- tools::makevars_user()
      makevars_size <<- file.size(makevars_during)
      stop("compilation failed")
    },
    rstan_options = function(...) FALSE
  )
  makevars_before <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  expect_error(cxxfunctionplus(), "compilation failed")
  # the user's Makevars was replaced by an empty file during compilation
  expect_length(makevars_during, 1)
  expect_equal(makevars_size, 0)
  expect_false(file.exists(makevars_during))
  # and R_MAKEVARS_USER is restored afterwards
  expect_identical(Sys.getenv("R_MAKEVARS_USER", unset = NA), makevars_before)
})
