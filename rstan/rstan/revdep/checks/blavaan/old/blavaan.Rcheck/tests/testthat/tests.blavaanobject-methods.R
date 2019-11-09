test_that("blavaan object methods work", {

  if(requireNamespace("rstan", quietly = TRUE) &
     requireNamespace("runjags", quietly = TRUE)){
    load(system.file("testdata", "sysdata.rda", package="blavaan"))

    # classes
    expect_equal(class(fitjags@external), "list")
    expect_equal(class(fitstan@external), "list")

    ## parameter summaries
    expect_equal(dim(parTable(fitjags)), c(10,20))
    expect_equal(dim(parTable(fitstan)), c(10,21))

    expect_equal(sum(fitjags@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitjags, 'psrf')))
    expect_equal(sum(fitstan@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitstan, 'psrf')))
    expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)
    expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstan)))
    
    expect_error(blavInspect(fitjags, 'blah'))

    ## fitMeasures
    expect_equal(length(fitMeasures(fitjags)),
                 length(fitMeasures(fitstan)))

    ## this is how summary() obtains its results, but have not figured out
    ## how to get S4 methods to directly work in testthat
    expect_equal(dim(parameterEstimates(fitjags)), c(10, 6))
    expect_equal(dim(parameterEstimates(fitstan)), c(10, 6))

    ## various blavInspect args
    expect_equal(length(blavInspect(fitjags, 'psrf')),
                 length(blavInspect(fitstan, 'psrf')))

    expect_equal(length(blavInspect(fitjags, 'neff')),
                 length(blavInspect(fitstan, 'neff')))

    expect_equal(length(blavInspect(fitjags, 'mcmc')),
                 length(blavInspect(fitstan, 'mcmc')))

    expect_equal(length(blavInspect(fitjags, 'start')),
                 length(blavInspect(fitstan, 'start')))

    expect_equal(dim(blavInspect(fitjags, 'hpd')),
                 dim(blavInspect(fitstan, 'hpd')))

    expect_equal(dim(standardizedposterior(fitjags)),
                 dim(standardizedposterior(fitstan)))
  }
    
})
