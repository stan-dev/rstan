x = getURL("http://192.168.0.3/stat141/Winter04/Homework/Solutions/index.html",
#                 userpwd = "name:passwd:",
#                  userpwd = "bob:",
                                     header = TRUE, verbose = TRUE,
                                     cookiefile = "/home/duncan/Rcookies",
                                     netrc = TRUE,
                                     maxredirs = as.integer(20),
                                     netrc.file = "/home2/duncan/.netrc1",
                                     followlocation = TRUE)

# This is no longer supported in libcurl.
# , passwdfunction = promptPassword)
