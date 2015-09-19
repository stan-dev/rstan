test_extract_sparse_parts <- function() {
  A <- rbind(
    c(19L, 27L,  0L,  0L),
    c( 0L,  0L,  0L,  0L),
    c( 0L,  0L,  0L, 52L),
    c(81L,  0L, 95L, 33L)
  )
  w <- c(19, 27, 52, 81, 95, 33)
  v <- c(1, 2, 4, 1, 3, 4)
  u <- c(1, 3, 3, 4, 7)
  parts <- rstan::extract_sparse_parts(A)
  checkEquals(parts$w, w)
  checkEquals(parts$v, v)
  checkEquals(parts$u, u)
  
  parts <- rstan::extract_sparse_parts(A * 1.0)
  checkEquals(parts$w, w)
  checkEquals(parts$v, v)
  checkEquals(parts$u, u)
  
  parts <- rstan::extract_sparse_parts(Matrix::Matrix(A))
  checkEquals(parts$w, w)
  checkEquals(parts$v, v)
  checkEquals(parts$u, u)
}

