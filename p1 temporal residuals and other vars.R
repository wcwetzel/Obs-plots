### plot1 temporal residuals and other vars
## 16 Nov 2011
library(bbmle)
library(ggplot2)

d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')

d1 = d1[!is.na(d1$galls2010),]
d1$mgalls = with(d1, (galls2011 + galls2010)/2 )
d1$vol = with(d1, d1/100 * d2/100 * h/100) # cubic meters
d1$density2011 = d1$galls2011/d1$vol


ggplot(data=d1, aes(x=galls2010, y=galls2011)) +
	geom_point() +
	stat_smooth()


m1 = mle2(galls2011 ~ dnbinom(mu = a + b * galls2010, size=k),
	start=list(a = 5, b = 0, k = 0.5), data=d1)
res1 = residuals(m1)

n = nrow(d1)
dists = pairdist(d1[,c('x','y')])
p1 = matrix(d1[,'galls2010'], nrow=n, ncol=n, byrow=FALSE)
p2 = matrix(d1[,'galls2010'], nrow=n, ncol=n, byrow=TRUE)
gdiff = abs(p1 - p2)

s1 = matrix(d1[,'vol'], nrow=n, ncol=n, byrow=FALSE)
s2 = matrix(d1[,'vol'], nrow=n, ncol=n, byrow=TRUE)
sdiff = abs(s1 - s2)


r1 = matrix(res1, nrow=n, ncol=n, byrow=FALSE)
r2 = matrix(res1, nrow=n, ncol=n, byrow=TRUE)
rdiff = abs(r1 - r2)


mantel(dists, rdiff)

mc.rdiff = mantel.correlog(D.eco=dists, D.geo=rdiff)
plot(mc.rdiff)


ggplot(data=d1, aes(x=meanmPa, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmeantot, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmeanclose, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmax, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vol, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=mgalls, y = res1)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmax, y = meanmPa)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmeantip, y = meanmPa)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

ggplot(data=d1, aes(x=vegmeantot, y = meanmPa)) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')


pca = with(d1[!is.na(d1$vegmax),], prcomp(~vegmax + vegmeanclose + 
	vegmeantip + meanmPa))
biplot(pca)

ggplot(data=d1, aes(y=res1[!is.na(vegmax)], x= pca$x[,1])) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')
	
	
ggplot(data=d1, aes(y=galls2011[!is.na(vegmax)], x= pca$x[,1])) +
	geom_point() +
	stat_smooth() +
	stat_smooth(method='lm')

mpc1 = mle2(galls2011[!is.na(vegmax)] ~ dnbinom())
