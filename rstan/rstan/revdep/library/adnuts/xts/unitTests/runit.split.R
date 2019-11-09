# 'f' is character, but length(f) > 1
test.split_character_f_not_endpoints <- function() {
  x <- .xts(1:5, 1:5)
  f <- letters[1:nrow(x)]
  checkIdentical(split(x,f), split(as.zoo(x),f))
}

