library(RCurl)

h = basicTextGatherer()
txt = getURL("http://192.168.0.3", header= TRUE, headerfunction = h[[1]], httpheader = c(Accept="text/html", Duncan="bob"), verbose = TRUE)
read.dcf(textConnection(paste(h$value(NULL)[-1], collapse="")))



