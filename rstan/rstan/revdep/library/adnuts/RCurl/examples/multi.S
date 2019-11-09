library(RCurl)
uris = c("http://www.omegahat.net/index.html", "http://www.omegahat.net/RecentActivities.html",
         "http://www.omegahat.net/RCurl/index.html", 
         "http://www.omegahat.net/RCurl/philosophy.xml",
         "http://www.omegahat.net/RCurl/philosophy.html")
write = multiTextGatherer(uris)
system.time(ans <- getURIAsynchronous(uris, write = write, perform = 1000))


system.time(invisible(getURI(uris, async = FALSE)))



atimes = sapply(1:40, function(i) system.time(getURI(uris)))
stimes = sapply(1:40, function(i) system.time(getURI(uris, async = FALSE)))

times = data.frame(user = c(atimes[1,], stimes[1,]),
                   system = c(atimes[2,], stimes[2,]),
                   elapsed = c(atimes[3,], stimes[3,]),
                   style = factor(c(rep("Asynchronous", 40), rep("Serial", 40))))

