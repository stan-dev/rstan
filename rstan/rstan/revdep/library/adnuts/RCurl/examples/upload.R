library(RCurl)
upFile = system.file("DESCRIPTION", package = "RCurl")
postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl",
         "fileData" = fileUpload(upFile), .opts = list(verbose = TRUE))


# Give the contents from R
postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl",
         "fileData" = fileUpload("", paste(readLines(upFile), collapse = "\n")),
         .opts = list(verbose = TRUE, header = TRUE))


postForm("http://eeyore.ucdavis.edu/cgi-bin/testForm1.pl",
         "fileData" = fileUpload("", paste(readLines(upFile), collapse = "\n"), "text/plain"),
         .opts = list(verbose = TRUE, header = TRUE))
