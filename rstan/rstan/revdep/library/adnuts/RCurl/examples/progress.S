library(RCurl)
getURL("http://www.omegahat.net/RCurl/index.html",  progressfunction = function(down, up) {print(down); print(up)}, noprogress = FALSE)
