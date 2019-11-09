### R code from vignette source 'xts.Rnw'

###################################################
### code chunk number 1: a
###################################################
require(xts)
data(sample_matrix)
class(sample_matrix)
str(sample_matrix)
matrix_xts <- as.xts(sample_matrix,dateFormat='Date')
str(matrix_xts)
df_xts <- as.xts(as.data.frame(sample_matrix),
 important='very important info!')
str(df_xts)


###################################################
### code chunk number 2: xtsconstructor
###################################################
xts(1:10, Sys.Date()+1:10)


###################################################
### code chunk number 3: xtsmethods (eval = FALSE)
###################################################
## matrix_xts['2007-03']


###################################################
### code chunk number 4: xtsmethods-hidden
###################################################
head(matrix_xts['2007-03'],5)
cat('...\n')


###################################################
### code chunk number 5: xtsmethods2 (eval = FALSE)
###################################################
## matrix_xts['/2007-01-07'] 


###################################################
### code chunk number 6: xtsmethods2-hidden
###################################################
matrix_xts['/2007-01-07']


###################################################
### code chunk number 7: xtsfirstandlast (eval = FALSE)
###################################################
## first(matrix_xts,'1 week')


###################################################
### code chunk number 8: xtsfirstandlast-hidden
###################################################
head(first(matrix_xts,'1 week'))


###################################################
### code chunk number 9: xtsfirstandlast2
###################################################
first(last(matrix_xts,'1 week'),'3 days')


###################################################
### code chunk number 10: indexClass
###################################################
indexClass(matrix_xts)
indexClass(convertIndex(matrix_xts,'POSIXct'))


###################################################
### code chunk number 11: xtsaxTicksByTime
###################################################
axTicksByTime(matrix_xts, ticks.on='months')


###################################################
### code chunk number 12: xtsplot
###################################################
plot(matrix_xts[,1],major.ticks='months',minor.ticks=FALSE,main=NULL,col=3)


###################################################
### code chunk number 13: asxtsreclass
###################################################
# using xts-style subsetting doesn't work on non-xts objects
sample_matrix['2007-06']
# convert to xts to use time-based subsetting
str(as.xts(sample_matrix)['2007-06'])

# reclass to get to original class back
str(reclass(as.xts(sample_matrix)['2007-06']))


###################################################
### code chunk number 14: usereclass
###################################################
z <- zoo(1:10,Sys.Date()+1:10)
# filter converts to a ts object - and loses the zoo class
(zf <- filter(z, 0.2))
class(zf)
# using Reclass, the zoo class is preserved
(zf <- Reclass(filter(z, 0.2)))
class(zf)


###################################################
### code chunk number 15: periodicity
###################################################
periodicity(matrix_xts)


###################################################
### code chunk number 16: endpoints
###################################################
endpoints(matrix_xts,on='months')
endpoints(matrix_xts,on='weeks')


###################################################
### code chunk number 17: toperiod
###################################################
to.period(matrix_xts,'months')
periodicity(to.period(matrix_xts,'months'))

# changing the index to something more appropriate
to.monthly(matrix_xts)


###################################################
### code chunk number 18: periodapply
###################################################
# the general function, internally calls sapply
period.apply(matrix_xts[,4],INDEX=endpoints(matrix_xts),FUN=max)


###################################################
### code chunk number 19: applymonthly
###################################################
# same result as above, just a monthly interface
apply.monthly(matrix_xts[,4],FUN=max)


###################################################
### code chunk number 20: periodsum
###################################################
# using one of the optimized functions - about 4x faster
period.max(matrix_xts[,4], endpoints(matrix_xts))


###################################################
### code chunk number 21: devtryxts
###################################################
period.apply


###################################################
### code chunk number 22: attributes
###################################################
str(attributes(matrix_xts))
str(xtsAttributes(matrix_xts))

# attach some attributes
xtsAttributes(matrix_xts) <- list(myattr="my meta comment")
attr(matrix_xts, 'another.item') <- "one more thing..."

str(attributes(matrix_xts))
str(xtsAttributes(matrix_xts))


###################################################
### code chunk number 23: subclass
###################################################
xtssubclass <- structure(matrix_xts, class=c('xts2','xts','zoo'))
class(xtssubclass)


