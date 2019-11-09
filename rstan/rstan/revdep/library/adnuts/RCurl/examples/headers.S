curl = getCurlHandle()
setCurlHeaders(c("Accept" = "text/html"), curl)
getURL("http://192.168.0.3", verbose = TRUE, curl = curl)
