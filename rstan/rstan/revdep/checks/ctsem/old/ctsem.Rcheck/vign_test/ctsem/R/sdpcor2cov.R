#' sdcor2cov
#' 
#' Converts a lower triangular matrix with standard deviations on the diagonal and partial correlations on
#' lower triangle, to a covariance (or cholesky decomposed covariance)
#' @param mat input square matrix with std dev on diagonal and lower tri of partial correlations.
#' @param cholesky Logical. To return the cholesky decomposition instead of full covariance, set to TRUE.
#' @examples
#' testmat <- diag(exp(rnorm(5,-3,2)),5) #generate arbitrary std deviations
#' testmat[row(testmat) > col(testmat)] <- runif((5^2-5)/2, -1, 1) 
#' print(testmat)
#' covmat <- sdpcor2cov(testmat) #convert to covariance
#' cov2cor(covmat) #convert covariance to correlation
#' @export
sdpcor2cov <- function(mat, cholesky=FALSE){ 

    ndim = ncol(mat);
    mcholcor=diag(0,ndim);
    mcholcor[1,1]=1;
    
    if(ndim > 1){
      for(coli in 1:ndim){
        for(rowi in coli:ndim){
          if(coli==1 && rowi > 1) mcholcor[rowi,coli] =  mat[rowi,coli]; 
          if(coli > 1){
            if(rowi == coli) mcholcor[rowi,coli] = prod(sqrt(1-mat[rowi,1:(coli-1)]^2));
            if(rowi > coli) mcholcor[rowi,coli] = mat[rowi,coli] * prod(sqrt(1-mat[rowi,1:(coli-1)]^2));
          }
        }
      }
    }
    mscale=diag(diag(mat))
    out= mscale %*% mcholcor
    if(!cholesky) out = out %*% t(out)
    return(out);
  }
 
