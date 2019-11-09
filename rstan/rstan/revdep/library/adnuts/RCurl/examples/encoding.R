library(RCurl)
Sys.setlocale(, "en_US.UTF-8")
x = getURL("http://www.omegahat.net/RCurl/index.html")
Encoding(x)

f = system.file("NAMESPACE", package = "RCurl")
x = postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl", "fileData" = fileUpload(f), .opts = list(header =TRUE))
x = postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl", "fileData" = fileUpload(f))
x = postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl", "fileData" = fileUpload(f),
              .opts = list(writefunction = function(x) {browser(); nchar(x)}))

# This determines the encoding from the HTTP header of the response
x = getURL("http://www.cl.cam.ac.uk/%7Emgk25/ucs/examples/UTF-8-demo.txt",
              header = TRUE,
               write = function(x) {print(Encoding(x))
                                    # print(nchar(x))  # this causes problems about invalid multibyte string. Why? Are we on a boundary
                                                       # i.e. the chunks come
                                    nchar(x, "bytes")})
