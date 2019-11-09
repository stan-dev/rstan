
doc = htmlTreeParse("http://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html", useInternal = TRUE)
types = xpathSApply(doc, "//h3", xmlValue)

types = gsub("^10\\.[0-9]+(\\.[0-9]+)? ", "", types)
types = types [ - grep("[0-9]xx$", types) ]

scan(con, what = c("integer", "character"))

status = gsub("^([0-9]+).*", "\\1", types)
name = gsub("[0-9]+ (.*)", "\\1", types)
name = gsub(" ", "_", name)

sQuote = function(x) paste("'", x, "'", sep = "")
cat("# Generated from inst/HTTPErrors/makeErrorClasses.R\n",
     "httpErrorClasses =\nc(", paste(sQuote( status ), sQuote(name), sep = " = ", collapse = ",\n   "), "\n)",
     file = "~/Projects/org/omegahat/R/RCurl/R/httpErrors.R")


