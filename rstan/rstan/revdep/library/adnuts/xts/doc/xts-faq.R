### R code from vignette source 'xts-faq.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("xts")
Sys.setenv(TZ="GMT")


###################################################
### code chunk number 2: xts-faq.Rnw:88-90 (eval = FALSE)
###################################################
## filenames <- c("a.csv", "b.csv", "c.csv")
## sample.xts <- as.xts(do.call("rbind", lapply(filenames, read.zoo)))


###################################################
### code chunk number 3: xts-faq.Rnw:106-107 (eval = FALSE)
###################################################
## lm(sample.xts[, "Res"] ~ sample.xts[, "ThisVar"] + sample.xts[, "ThatVar"])


###################################################
### code chunk number 4: xts-faq.Rnw:110-111 (eval = FALSE)
###################################################
## with(sample.xts, lm(Res ~ ThisVar + ThatVar))


###################################################
### code chunk number 5: xts-faq.Rnw:118-121
###################################################
sample.xts <- xts(c(1:3, 0, 0, 0), as.POSIXct("1970-01-01")+0:5)
sample.xts[sample.xts==0] <- NA
cbind(orig=sample.xts, locf=na.locf(sample.xts))


###################################################
### code chunk number 6: xts-faq.Rnw:128-130
###################################################
data(sample_matrix)
sample.xts <- xts(1:10, seq(as.POSIXct("1970-01-01"), by=0.1, length=10))


###################################################
### code chunk number 7: xts-faq.Rnw:138-140
###################################################
options(digits.secs=3)
head(sample.xts)


###################################################
### code chunk number 8: xts-faq.Rnw:151-154
###################################################
dt <- as.POSIXct("2012-03-20 09:02:50.001")
print(as.numeric(dt), digits=20)
sprintf("%20.10f", dt)


###################################################
### code chunk number 9: xts-faq.Rnw:163-164 (eval = FALSE)
###################################################
## sample.xts.2 <- xts(t(apply(sample.xts, 1, myfun)), index(sample.xts))


###################################################
### code chunk number 10: xts-faq.Rnw:172-177
###################################################
sample.xts <- xts(1:50, seq(as.POSIXct("1970-01-01"),
  as.POSIXct("1970-01-03")-1, length=50))
apply.daily(sample.xts, mean)
period.apply(sample.xts, endpoints(sample.xts, "days"), mean)
period.apply(sample.xts, endpoints(sample.xts, "hours", 6), mean)


###################################################
### code chunk number 11: xts-faq.Rnw:185-186 (eval = FALSE)
###################################################
## apply.daily(sample.xts['T06:00/T17:00',], mean)


###################################################
### code chunk number 12: xts-faq.Rnw:196-209
###################################################
sample.xts <- xts(1:6, as.POSIXct(c("2009-09-22 07:43:30",
  "2009-10-01 03:50:30", "2009-10-01 08:45:00", "2009-10-01 09:48:15",
  "2009-11-11 10:30:30", "2009-11-11 11:12:45")))
# align index into regular (e.g. 3-hour) blocks
aligned.xts <- align.time(sample.xts, n=60*60*3)
# apply your function to each block
count <- period.apply(aligned.xts, endpoints(aligned.xts, "hours", 3), length)
# create an empty xts object with the desired regular index
empty.xts <- xts(, seq(start(aligned.xts), end(aligned.xts), by="3 hours"))
# merge the counts with the empty object
head(out1 <- merge(empty.xts, count))
# or fill with zeros
head(out2 <- merge(empty.xts, count, fill=0))


###################################################
### code chunk number 13: xts-faq.Rnw:219-220 (eval = FALSE)
###################################################
## sample.xts <- as.xts(transform(sample.xts, ABC=1))


###################################################
### code chunk number 14: xts-faq.Rnw:224-225 (eval = FALSE)
###################################################
## indexTZ(sample.xts) <- Sys.getenv("TZ")


###################################################
### code chunk number 15: xts-faq.Rnw:237-238 (eval = FALSE)
###################################################
## sample.xts[sample.xts$Symbol == "AAPL" & index(sample.xts) == as.POSIXct("2011-09-21"),]


###################################################
### code chunk number 16: xts-faq.Rnw:241-242 (eval = FALSE)
###################################################
## sample.xts[sample.xts$Symbol == "AAPL"]['2011-09-21']


###################################################
### code chunk number 17: xts-faq.Rnw:249-253
###################################################
data(sample_matrix)
sample.xts <- as.xts(sample_matrix)
wday.xts <- sample.xts[.indexwday(sample.xts) %in% 1:5]
head(wday.xts)


###################################################
### code chunk number 18: xts-faq.Rnw:266-268
###################################################
Data <- data.frame(timestamp=as.Date("1970-01-01"), obs=21)
sample.xts <- xts(Data[,-1], order.by=Data[,1])


###################################################
### code chunk number 19: xts-faq.Rnw:272-275
###################################################
Data <- data.frame(obs=21, timestamp=as.Date("1970-01-01"))
sample.xts <- xts(Data[,!grepl("timestamp",colnames(Data))],
  order.by=Data$timestamp)


###################################################
### code chunk number 20: xts-faq.Rnw:288-291 (eval = FALSE)
###################################################
## x1 <- align.time(xts(Data1$obs, Data1$timestamp), n=600)
## x2 <- align.time(xts(Data2$obs, Data2$timestamp), n=600)
## merge(x1, x2)


###################################################
### code chunk number 21: xts-faq.Rnw:295-301
###################################################
data(sample_matrix)
sample.xts <- as.xts(sample_matrix)
sample.xts["2007-01"]$Close <- sample.xts["2007-01"]$Close + 1
#Warning message:
#In NextMethod(.Generic) :
#  number of items to replace is not a multiple of replacement length


###################################################
### code chunk number 22: xts-faq.Rnw:314-315 (eval = FALSE)
###################################################
## sample.xts["2007-01",]$Close <- sample.xts["2007-01"]$Close + 1


###################################################
### code chunk number 23: xts-faq.Rnw:323-324 (eval = FALSE)
###################################################
## sample.xts["2007-01","Close"] <- sample.xts["2007-01","Close"] + 1


