# this wrong model is supposed to give very long
# error message (more than 1000 characters)
parameters {
  cov_matrix[3] y;
} 
model {
  y ~ normal(0, 1); 
}
