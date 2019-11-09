curl = getCurlHandle()
x = getURL("http://192.168.0.3", verbose = TRUE, curl = curl)
getCurlInfo(.Call("R_getCURLInfoEnum"), curl)


