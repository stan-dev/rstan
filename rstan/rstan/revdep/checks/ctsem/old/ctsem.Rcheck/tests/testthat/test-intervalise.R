library(ctsem)
library(testthat)

context("intervalise")

test_that("intervals", {
  data('longexample')
  
  #Then convert to wide format
  wideexample <- ctLongToWide(datalong = longexample, id = "id",
                              time = "time", manifestNames = c("Y1", "Y2", "Y3"),
                              TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
  
  #Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
  wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3,
                        n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"),
                        TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )
  
  dt <- matrix(c(.89, .53, .001, .210), nrow=2,ncol=2,
            ,dimnames = list(c("1", "2"), c("dT1", "dT2")))
  expect_equal(wide[,c('dT1','dT2')], dt)
})
