## partial mantel with distances, galls, and plant size
# 10 Nov 2011

library(spatstat)
library(vegan)
library(compiler)


#d1 = read.csv(
#	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
#d2 = d1

d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol




n = nrow(d2)
dists = pairdist(d2[,c('x','y')])
p1 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=FALSE)
p2 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=TRUE)
gdiff = abs(p1 - p2)

s1 = matrix(d2[,'vol'], nrow=n, ncol=n, byrow=FALSE)
s2 = matrix(d2[,'vol'], nrow=n, ncol=n, byrow=TRUE)
sdiff = abs(s1 - s2)

mantel.partial(dists, gdiff, sdiff)


