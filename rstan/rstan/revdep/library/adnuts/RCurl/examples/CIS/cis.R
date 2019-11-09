# From Roger Koenker.
"http://www.statindex.org/CIS/CIS/CIS/psqlQuery/Searcher?>>> > skip=+&authed=>>> > +&authed=a&authed=e&authorstring=koenker&titlestring=rank&keywordstrin>>> > g=&jnt=jp&jnamestring=&pt=+&pt=b&pt=j&pt=p&pt=+&yearbeg=&yearend=&startDisplay=1&endDisplay=50&fmt=extended

scis <- "http://www.statindex.org/CIS/CIS/CIS/psqlQuery/Searcher?skip=+&authed=+&authed=a&authed=e&authorstring=koenker&titlestring=rank&keywordstring=&jnt=jp&jnamestring=&pt=+&pt=b&pt=j&pt=p&pt=+&yearbeg=&yearend=&startDisplay=1&endDisplay=50&fmt=extended" 



h = getCurlHandle(cookiejar = "/tmp/MyCookies")

# XXX Fill in user id and password here.
postForm("https://www.amstat.org/membersonly/index.cfm?fuseaction=login", txtUser = "106823", txtPassword = "Jasper8Hazel", curl = h)

#getURL("https://www.amstat.org/membersonly/index.cfm?fuseaction=CISWeb")

# This gets us back the _ZopeId cookie
v = getURI("https://www.statindex.org/psqlQuery", curl = h)

u = "http://www.statindex.org/psqlQuery"
u = "http://www.statindex.org/CIS/psqlQuery/Searcher"
u = "http://www.statindex.org/psqlQuery/Searcher"


newHandle = getCurlHandle(cookie = '_ZopeId="44748384A25VQtBtVWY"')

cookie = '_ZopeId="59407264A25VRGMoY0Q"'

tmp = postForm(u,
          skip = "+", authed = "+", authed = "a", authorstring = "Ihaka", titlestring = "",
          keywordstring ="", jnt = "jp", jnamestring = "", pt= "+", pt = "b", pt = "j", pt = "p", pt = "+",
          yearbeg = "", yearend = "",
          startDisplay = "1", endDisplay = "50", fmt = "extended",
          .opts = list(verbose = TRUE), curl = h)




# cookie = '_ZopeId="95878868A25Ttl.5Nd0"'

getForm("https://www.statindex.org/psqlQuery")


###########################################################################


