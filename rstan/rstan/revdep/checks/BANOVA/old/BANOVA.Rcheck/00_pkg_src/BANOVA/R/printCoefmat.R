printCoefmat <-
function (coef_table){
  print(data.frame(coef_table))
  cat('---')
  cat('\nSignif. codes: 0 *** 0.001 ** 0.01 * 0.05 . 0.1 1 \n')
}
