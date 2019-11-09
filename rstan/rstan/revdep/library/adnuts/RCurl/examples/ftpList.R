
library(RCurl)
url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/"
fileNames = getURL(url, .opts = list(customrequest = "NLST *.gz") )
fileNames = strsplit(fileNames, "\\\r\\\n")[[1]]

# Now you can download these directly but you have to deal
# with the compression. This is possible with RCurl by
# specifying a binary reader.

# Alternatively, get the entiore directory listing and have to pull
# the names out of this much bigger download.
z = getURL(url, .opts = list(customrequest = "LIST *.gz") )
