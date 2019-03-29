get_stream <- function() {
  message("The 'get_stream' function is deprecated and does not do anything useful.",
          "The output stream is now handled by Rcpp.")
  return(.Call('get_stream_'))
}



