# prelim look at 2010 p1 galls
# 7 Apr 2012


### some early results
# ptism rate doesn't vary with gall abundance, or maybe has a weak
	# humped shaped relationship, with a peak ptism at intermediate
	# gall densities.
# predation rate seems to increase slightly with gall abundance
# survival rate doesn't change or has a very slight decrease

library(ggplot2)
library(bbmle)


data = read.csv('~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis.csv')
data = data[order(data$plant),] # sort data on plant ID number

# tabulate gall outcomes by plant
gall.count = table(data$plant)
healthy.count = as.vector(by(data$HEALTHY, data$plant, sum))
ptoid.count = as.vector(by(data$PARASITIZED, data$plant, sum))
pred.count = as.vector(by(data$PREDATION, data$plant, sum))
unk.count = as.vector(by(data$UNK_MORT, data$plant, sum))

# make new data frame in which each row is a plant
data2 = data.frame(
	plant = dimnames(gall.count)[[1]], 
	galls = as.vector(gall.count),
	healthy = healthy.count,
	ptism = ptoid.count,
	pred = pred.count,
	unk = unk.count
	)

# add per plant rates to data frame
data2$s = with(data2, healthy/galls)
data2$ptism.rate = with(data2, ptism/galls)
data2$pred.rate = with(data2, pred/galls)
data2$unk.rate = with(data2, unk/galls)

# survival excluding ptism and predation
# is there DD on plant ignoring predation/parasitism?
data2$s2 = with(data2,  healthy / (healthy + unk))
#data2$s2[is.nan(data2$s2)] = 0 # don't want to do this. should be NaN


#--------------------- plotting ----------------------------#

hist(data2$galls)
hist(data2$healthy)
hist(data2$ptism)
hist(data2$pred)
hist(data2$unk)

# survival rate as a func of # galls
ggplot(aes(x=galls, y=s), data=data2) +
	geom_point() +
	stat_smooth()
	
# prop of healthy to healthy+unk
ggplot(aes(x=galls, y=s2), data=data2) +
	geom_point() +
	stat_smooth()

# ptism rate as a func of # galls
ggplot(aes(x=galls, y=ptism.rate), data=data2) +
	geom_point() +
	stat_smooth()

# pred rate as a func of # galls
ggplot(aes(x=galls, y=pred.rate), data=data2) +
	geom_point() +
	stat_smooth()


# log scale #
# survival rate as a func of # galls
ggplot(aes(x=log(galls), y=s), data=data2) +
	geom_point() +
	stat_smooth()
	
# prop of healthy to healthy+unk
ggplot(aes(x=log(galls), y=s2), data=data2) +
	geom_point() +
	stat_smooth()

# ptism rate as a func of # galls
ggplot(aes(x=log(galls), y=ptism.rate), data=data2) +
	geom_point() +
	stat_smooth()

# pred rate as a func of # galls
ggplot(aes(x=log(galls), y=pred.rate), data=data2) +
	geom_point() +
	stat_smooth()



#--------------- some models ---------------#

## survival ##
m0.healthy = mle2(healthy ~ dbinom(size=galls, prob=plogis(p)), data=data2, 
	start=list(p=0))

m1.healthy = mle2(healthy ~ dbinom(size=galls, prob=plogis(p + beta * galls)), 
	data=data2, start=list(p=0, beta=0))

m2.healthy = mle2(healthy ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2)), data=data2, start=list(p=0, beta=0, beta2=0))

m3.healthy = mle2(healthy ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2 + beta3 * galls^3)), data=data2, 
	start=list(p=0.1, beta=0.089, beta2=-0.057, beta3=0.002486))

AICtab(m0.healthy, m1.healthy, m2.healthy, m3.healthy)

predicted.p0 = predict(m0.healthy)/data2$galls
predicted.p1 = predict(m1.healthy)/data2$galls
predicted.p2 = predict(m2.healthy)/data2$galls
predicted.p3 = predict(m3.healthy)/data2$galls


plot(s ~ galls, data=data2, pch=20)
points(predicted.p0 ~ galls, data=data2, pch='0', col='red')
points(predicted.p1 ~ galls, data=data2, pch='1', col='blue')
points(predicted.p2 ~ galls, data=data2, pch='2', col='green')
points(predicted.p3 ~ galls, data=data2, pch='3', col='yellow')





## parasitism ##
m0.ptism = mle2(ptism ~ dbinom(size=galls, prob=plogis(p)), data=data2, 
	start=list(p=0))

m1.ptism = mle2(ptism ~ dbinom(size=galls, prob=plogis(p + beta * galls)), 
	data=data2, start=list(p=0, beta=0))

m2.ptism = mle2(ptism ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2)), data=data2, start=list(p=0, beta=0, beta2=0))

m3.ptism = mle2(ptism ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2 + beta3 * galls^3)), data=data2, 
	start=list(p=0.1, beta=0.089, beta2=-0.057, beta3=0.002486))

AICtab(m0.ptism, m1.ptism, m2.ptism, m3.ptism)

predicted.p = predict(m2.ptism)/data2$galls

plot(ptism.rate ~ galls, data=data2)
points(predicted.p ~ galls, data=data2, pch='G')

## predation ##
m0.pred = mle2(pred ~ dbinom(size=galls, prob=plogis(p)), data=data2, 
	start=list(p=0))

m1.pred = mle2(pred ~ dbinom(size=galls, prob=plogis(p + beta * galls)), 
	data=data2, start=list(p=0, beta=0))

m2.pred = mle2(pred ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2)), data=data2, start=list(p=0, beta=0, beta2=0))

m3.pred = mle2(pred ~ dbinom(size=galls, prob=plogis(p + beta * galls +
	beta2 * galls^2 + beta3 * galls^3)), data=data2, 
	start=list(p=0.1, beta=0.089, beta2=-0.057, beta3=0.002486))

AICtab(m0.pred, m1.pred, m2.pred, m3.pred)

predicted.p0 = predict(m0.pred)/data2$galls
predicted.p1 = predict(m1.pred)/data2$galls
predicted.p2 = predict(m2.pred)/data2$galls
predicted.p3 = predict(m3.pred)/data2$galls


plot(pred.rate ~ galls, data=data2, pch=20)
points(predicted.p0 ~ galls, data=data2, pch='0', col='red')
points(predicted.p1 ~ galls, data=data2, pch='1', col='blue')
points(predicted.p2 ~ galls, data=data2, pch='2', col='green')
points(predicted.p3 ~ galls, data=data2, pch='3', col='yellow')


