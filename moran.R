### moran's I for residuals of temporal
# 10 Nov 2011

library(spatstat)
library(vegan)
library(compiler)
library(bbmle)
library(ape)
library(spdep)




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


m2 = mle2(galls2011 ~ dnbinom(mu = a + b * galls2010, size=k),
	start=list(a = 5, b = 0, k = 0.5), data=d2)
res2 = residuals(m2)

m1 = mle2(galls2011 ~ dnbinom(mu = a + b * galls2010, size=k),
	start=list(a = 5, b = 0, k = 1), data=d1)
res1 = residuals(m1)

plot(d2$y ~ d2$x, cex=res-min(res) + 0.25)
points(-0.7,-0.9, cex = 0-min(res) + 0.25, pch=1)
plot(d1$y ~ d1$x, cex=res-min(res) + 0.25)
points(-0.7,-0.9, cex = 0-min(res) + 0.25, pch=1)

# for plot2
# first create a matrix of inverse distance weights
dists = as.matrix(dist(cbind(d2$x, d2$y)))
dists.inv2 = 1/dists
diag(dists.inv2) = 0
dists.inv2[1:5,1:5]
Moran.I(d2$galls2010, dists.inv2)
Moran.I(d2$galls2011, dists.inv2)

# can use Moran's I with binned distances
dists.bin1 = dists>0 & dists <= 1
Moran.I(d2$galls2010, dists.bin1)
Moran.I(d2$galls2011, dists.bin1)
dists.bin2 = dists>0 & dists <= 2
Moran.I(d2$galls2010, dists.bin2)
Moran.I(d2$galls2011, dists.bin2)
dists.bin3 = dists>0 & dists <= 3
Moran.I(d2$galls2010, dists.bin3)
Moran.I(d2$galls2011, dists.bin3)

# Moran's I on residuals
Moran.I(res2, dists.inv2)

# for plot1
# first create a matrix of inverse distance weights
dists = as.matrix(dist(cbind(d1$x, d1$y)))
dists.inv1 = 1/dists
diag(dists.inv1) = 0
dists.inv1[1:5,1:5]
Moran.I(d1$galls2010, dists.inv1)
Moran.I(d1$galls2011, dists.inv1)

# can use Moran's I with binned distances
dists.bin1 = dists>0 & dists <= 1
Moran.I(d1$galls2010, dists.bin1)
Moran.I(d1$galls2011, dists.bin1)
dists.bin2 = dists>0 & dists <= 2
Moran.I(d1$galls2010, dists.bin2)
Moran.I(d1$galls2011, dists.bin2)
dists.bin3 = dists>0 & dists <= 3
Moran.I(d1$galls2010, dists.bin3)
Moran.I(d1$galls2011, dists.bin3)

# Moran's I on residuals
Moran.I(res1, dists.inv1)

