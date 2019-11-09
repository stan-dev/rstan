rowMatch <-
function(vector, matrix){
  if (length(vector) == 1) return(match(vector, matrix))
  for (i in 1:nrow(matrix)){
    if (sum(vector == matrix[i,]) == length(vector))
      return(i)
  }
  return(NA)
}
