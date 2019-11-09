# Winnie on home machine.
getURI("https://192.168.0.2/stat141/Winter04/Homework/Solutions/",
       userpwd = "bob:ComeIn",
       ssl.verifyhost = as.integer(0), ssl.verifypeer = as.integer(0))


# Omitting the verifyhost causes
# SSL: certificate subject name 'localhost.localdomain' does not match target host name '192.168.0.2'

# However, specifying 1 as its value works even though it should be verifying the host.
#  2 is the default and forces the error.

