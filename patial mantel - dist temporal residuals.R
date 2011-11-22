### partial mantel for distances, residuals of 2010-1
# 10 Nov 2011

library(spatstat)
library(vegan)
library(compiler)
library(bbmle)

d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
d1 = d1[!is.na(d1$galls2010),]
d1$mgalls = with(d1, (galls2011 + galls2010)/2 )
d1$vol = with(d1, d1/100 * d2/100 * h/100) # cubic meters
d1$density = d1$galls2011/d1$vol


d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol

#d2=d1

m = mle2(galls2011 ~ dnbinom(mu = a + b * galls2010, size=k),
	start=list(a = 5, b = 0, k = 0.5), data=d2)

res = residuals(m)


n = nrow(d2)
dists = pairdist(d2[,c('x','y')])
p1 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=FALSE)
p2 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=TRUE)
gdiff = abs(p1 - p2)

s1 = matrix(d2[,'vol'], nrow=n, ncol=n, byrow=FALSE)
s2 = matrix(d2[,'vol'], nrow=n, ncol=n, byrow=TRUE)
sdiff = abs(s1 - s2)


r1 = matrix(res, nrow=n, ncol=n, byrow=FALSE)
r2 = matrix(res, nrow=n, ncol=n, byrow=TRUE)
rdiff = abs(r1 - r2)
rprod = r1 * r2


mantel(dists, rdiff)

mc.rdiff = mantel.correlog(D.eco=dists, D.geo=rdiff)
plot(mc.rdiff)
#title(main = 'Plot 2: Mantel correlogram')
