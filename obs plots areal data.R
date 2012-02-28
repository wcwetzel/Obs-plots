### Analyzing obs plots with plants as areal data
## 20 Feb 2012

# load packages
library(spdep)
library(bbmle)
library(sp)
library(pgirmess)

# load data
d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol
d2$change = d2$galls2011 - d2$galls2010

# turn x-y coordinates into a square distance matrix
dists = as.matrix(dist(cbind(d2$x,d2$y), diag=TRUE, upper=TRUE))

# now convert distance matrix into weights matrix
# with negative exp decline in weight with distance
curve(exp(-x * 1), 0, max(dists))
abline(h=0.05)
abline(h=0.1)
abline(h=0.25)
nexpW = exp(-dists)

# mat2listw converts a square spatial weights matrix into 
# a spdep weights list object
nelistW = mat2listw(nexpW, row.names=d2$tag)

# now try a neg exp weight with a cut off at 0.15
# with negative exp decline in weight with distance
nexpcoW = nexpW
nexpcoW[nexpcoW < 0.25] = 0

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

moran.mc(uncorr, necolistW, nsim=4999)
moran.mc(d2$galls2011, necolistW, nsim=4999)
moran.mc(d2$galls2010, necolistW, nsim=4999)
moran.mc(d2$mgalls, necolistW, nsim=4999)
mcm.neco = moran.mc(d2$change, necolistW, nsim=4999)
mcm.ne = moran.mc(d2$change, nelistW, nsim=4999)



hikingboots = mcm.ne$res[-length(mcm.ne$res)]
bootdens = density(hikingboots)
hist(hikingboots, prob=TRUE)
lines(bootdens)
abline(v=mcm.ne$statistic,lwd=2,col=4)



# moran scatter plot:
# this shows how similar a point is to its neighbors
# x-axis is d2$change
# y-axis is a weighted average of neighbors
moran.plot(d2$change, nelistW)
# there seems to be one far out influential point
# but even with that point removed, there is a strong trend
# for points with a bigger change
# to have neighbors with a bigger change


# what about using residuals instead of raw change?

model = glm(galls2011 ~ galls2010, family=poisson, data=d2)
res = resid(model)
mcm.ne.res = moran.mc(res, nelistW, nsim=4999)
moran.plot(res, nelistW)


hikingboots2 = mcm.ne.res$res[-length(mcm.ne.res$res)]
bootdens2 = density(hikingboots2)
hist(hikingboots2, prob=TRUE)
lines(bootdens2)
abline(v=mcm.ne.res$statistic,lwd=2,col=4)

## look for non-stationarity ##
plot(galls2011 ~ x, data=d2)
plot(galls2010 ~ x, data=d2)
plot(change ~ x, data=d2)

plot(res ~ x, data=d2)
summary(lm(res~ d2$x))
abline(lm(res~ d2$x))

plot(res ~ y, data=d2)
summary(lm(res~ d2$y))
abline(lm(res~ d2$y))
# It looks as though change in galls increases with x coord!

# residuals from a model that has x coords in it too
model2 = glm(galls2011 ~ galls2010 + x, family=poisson, data=d2)
res2 = resid(model2)
mcm.ne.res2 = moran.mc(res2, nelistW, nsim=4999)
moran.plot(res2, nelistW)



# spatial correlograms
#neighbors = listw2sn(nelistW)
neighbors = dnearneigh(coordinates(d2), d1=0, d2=1.5, row.names=d2$tag)
summary(neighbors)
plot(neighbors, coordinates(d2))
cor1 = sp.correlogram(neighbors, var=d2$change, order=8)
print(cor1, p.adj.method='holm')
plot(cor1)


# pgirmess correlog
cor2 = correlog(coordinates(d2), d2$change)
plot(cor2)

cor3 = correlog(coordinates(d2), d2$galls2011)
plot(cor3)


