library(RCurl)

# From Seth Falcon.

getDownloadSize = function (url)
{
  h <- basicTextGatherer()
  junk <- getURI(url, headerfunction = h$update, header = TRUE, nobody = TRUE)
  h <- h$value()
  ##parseContentLength(h)
  h
}

u = "http://cran.fhcrc.org/src/contrib/PACKAGES.gz"

# If we do 100 iterations, we fail on the 96th with an error
# about not being able to resolve cran.fchrc.org.
# This is on my linux box at home.  If we sleep for a second between
# each call, all is well. If we do just 90 even without sleeping, all is well.
# So looks like there is some maximum number of requests per time period on that
# DNS perhaps for the same machine.....
for (i in 1:90) {
   print(i)
#   Sys.sleep(1)
   jj = getDownloadSize(u)
}

