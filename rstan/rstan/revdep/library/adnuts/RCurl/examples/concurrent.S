# See concurrent.xml for more details

library(RCurl)

getURIs =
function(uris, ..., multiHandle = getCurlMultiHandle(), .perform = TRUE)
{
  content = list()
  curls = list()

  for(i in uris) {
    curl = getCurlHandle()
    content[[i]] = basicTextGatherer()
    opts = curlOptions(URL = i, writefunction = content[[i]]$update, ...)    
    curlSetOpt(.opts = opts, curl = curl)
    multiHandle = push(multiHandle, curl)
  }

  
  if(.perform) {
     complete(multiHandle)
     lapply(content, function(x) x$value())
   } else {
     return(list(multiHandle = multiHandle, content = content))
   }
}


uris = c("http://www.omegahat.net/index.html", 
         "http://www.r-project.org/src/contrib/PACKAGES.html",
         "http://developer.r-project.org/index.html", 
         "http://www.slashdot.org/philosophy.xml",
         "http://fxfeeds.mozilla.org/rss20.xml",
         "http://www.nytimes.com/index.html")



#atimes = sapply(1:40, function(i) system.time(getURIs(uris, maxconnects = 100)))
#stimes = sapply(1:40, function(i) system.time(lapply(uris, getURI)))

n = 100
  # do one function call for each URI
serialTimes = replicate(n, system.time(lapply(uris, getURI)))

  # do all in function call but sequentially
vectorizedTimes = replicate(n, system.time(getURI(uris, async = FALSE)))

  # use the asynchronous approach
performingAsyncTimes = replicate(n, 
                   system.time({
                     getURIs(uris)
                   }))

   # And do the asynchronous approach more manually
asyncTimes = replicate(n, 
                   system.time({
                     z = getURIs(uris, .perform = FALSE)
                     complete(z$multiHandle)
                     lapply(z$content, function(x) x$value())
# Can fetch the download times of the individual documents                     
#sapply(z$multiHandle@subhandles, function(x) getCurlInfo(x)$total.time)                     
                   }))



timings =
  data.frame(user = c(performingAsyncTimes[1,], vectorizedTimes[1,], serialTimes[1,]),
           system = c(performingAsyncTimes[2,], vectorizedTimes[2,], serialTimes[2,]),
           elapsed = c(performingAsyncTimes[3,], vectorizedTimes[3,], serialTimes[3,]),
           mode = factor(rep(c("Asynchronous", "Vectorized", "Serial"), rep(n, 3))))


