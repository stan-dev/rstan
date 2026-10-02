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
