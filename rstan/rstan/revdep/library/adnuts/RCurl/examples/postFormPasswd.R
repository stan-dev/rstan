# Perl script (and HTML form for testing in the browser) taken from
#   http://www.elated.com/articles/form-validation-with-perl-and-cgi/


# Provide the login & password directly
postForm("http://www.omegahat.net/RCurl/testPassword/form_validation.cgi",
           your_name = "Duncan",
           your_age = "35-55",
           your_sex = "m",
           submit = "submit",
           .opts = list(userpwd = "bob:welcome"))

# Get the login & password in ~/.netrc
postForm("http://www.omegahat.net/RCurl/testPassword/form_validation.cgi",
           your_name = "Duncan",
           your_age = "35-55",
           your_sex = "m",
           submit = "submit",
          .opts = list(netrc = TRUE))

# Get the login & password from a different netrc file

postForm("http://www.omegahat.net/RCurl/testPassword/form_validation.cgi",
           your_name = "Duncan",
           your_age = "35-55",
           your_sex = "m",
           submit = "submit",
           .opts = list(netrc = TRUE,
                        netrc.file = "/Users/duncan/Projects/org/omegahat/R/RCurl/inst/examples/omg.netrc"))
