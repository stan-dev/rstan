summer =
function()
{
  total = 0.0
  numValues = 0

  list(read = function(txt) {
                 con = textConnection(txt)
                 on.exit(close(con))
                 els = scan(con)
                 numValues <<- numValues + length(els)
                 total <<- total + sum(els)

                 ""
              },
       result = function() c(total = total, numValues = numValues))
}

s = summer()

getURL("http://www.omegahat.net/RCurl/matrix.data", write = chunkToLineReader(s$read)$read)
