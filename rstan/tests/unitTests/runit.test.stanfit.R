.setUp <- function() {
  csv_fname <- 'teststanfit.csv'
  model_code <- 
    "parameters { 
       real y[2];
     } 
     transformed parameters {
       real y2[2, 2];
       y2[1, 1] <- y[1]; 
       y2[1, 2] <- -y[1]; 
       y2[2, 1] <- -y[2]; 
       y2[2, 2] <- y[2]; 
     } 
     model {
       y ~ normal(0, 1);
     } 
     generated quantities {
       real y3[3, 2, 3];
       y3[1,1,1] <- 1;   y3[1,1,2] <- 2;   y3[1,1,3] <- 3;  
       y3[1,2,1] <- 4;   y3[1,2,2] <- 5;   y3[1,2,3] <- 6;
       y3[2,1,1] <- 7;   y3[2,1,2] <- 8;   y3[2,1,3] <- 9;  
       y3[2,2,1] <- 10;  y3[2,2,2] <- 11;  y3[2,2,3] <- 12;
       y3[3,1,1] <- 13;  y3[3,1,2] <- 14;  y3[3,1,3] <- 15;  
       y3[3,2,1] <- 16;  y3[3,2,2] <- 17;  y3[3,2,3] <- 18;
     } 
    "
    fit <- stan(model_code = model_code, 
                iter = 100, 
                chains = 1,
                sample_file = csv_fname)
} 
 

test_output_csv <- function() {
  csv_fname <- 'teststanfit.csv'
  checkTrue(file.exists(csv_fname))
  d <- read.csv(file = csv_fname, comment.char = '#',
                header = TRUE)
  y_names <- paste0("y.", 1:2)
  cat2 <- function(a, b) paste0(a, ".", b)
  y2_names <- paste0("y2.", t(outer(1:2, 1:2, cat2)))
  y3_names <- paste0("y3.", t(outer(as.vector(t(outer(1:3, 1:2, cat2))), 1:3, cat2)))

  # FIXME, uncomment the follwing line
  # checkEquals(colnames(d), c("lp__", "treedepth__", "stepsize__", y_names, y2_names, y3_names))
  iter1 <- unlist(d[1,])
  # to check if the order is also row-major
  checkEquals(iter1["y.1"], iter1["y2.1.1"], checkNames = FALSE)
  checkEquals(iter1["y.1"], -iter1["y2.1.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], iter1["y2.2.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], -iter1["y2.2.1"], checkNames = FALSE)
  checkEquals(iter1[y3_names], 1:18, checkNames = FALSE)
} 

.tearDown <- function() {
  csv_fname <- 'teststanfit.csv'
  unlink(csv_fname)
}

