### Analyzing obs plots with plants as areal data
## 20 Feb 2012

# load packages
library(spdep)
library(bbmle)
library(sp)

# load data
d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol

# turn x-y coordinates into a square distance matrix
dists = as.matrix(dist(cbind(d2$x,d2$y), diag=TRUE, upper=TRUE))

# now convert distance matrix into weights matrix
# with negative exp decline in weight with distance
curve(exp(-x * 1), 0, max(dists))
nexpW = exp(-dists)

# mat2listw converts a square spatial weights matrix into 
# a spdep weights list object
nelistW = mat2listw(nexpW, row.names=d2$tag)

# now try a neg exp weight with a cut off at 0.05
# with negative exp decline in weight with distance
nexpcoW = nexpW
nexpcoW[nexpcoW < 0.05] = 0

# mat2listw converts a square spatial weights matrix into 
# a spdep weights list object
necolistW = mat2listw(nexpcoW, row.names=d2$tag)




# some things we can do with spdep spatial weights objects
summary(nelistW)
summary(necolistW)

coordinates(d2) = cbind(d2$x, d2$y)

plot(nelistW, coordinates(d2))
plot(necolistW, coordinates(d2))


# OK now Moran's I
uncorr = rpois(nrow(d2), lambda=mean(d2$mgalls))
mt.ne = moran.test(d2$galls2011, listw=nelistW)
mt.neco = moran.test(d2$galls2011, listw=necolistW)
moran.test(d2$mgalls, listw=necolistW)
moran.test(d2$galls2010, listw=necolistW)
moran.test(uncorr, listw=necolistW)
# it's not working. i'm getting very sig autocorr for
# a random distribution of galls. crazy!





