# R function to call your compiled code

checkOrder <- function(x) {
  .Call('check_order', x, TRUE, TRUE)
}
