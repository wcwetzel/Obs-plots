## gall fate, size, and flowers ##
# 6 May 2012

library(ggplot2)
library(bbmle)

data2 = read.csv('~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis flws.csv')
# oops this doesn't have mean gall size

# let's put it in
data = read.csv('~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis.csv')
data = data[order(data$plant),] # sort data on plant ID number


# mean gall size by plant
diam.mean = as.vector(by(data$gall_diameter, data$plant, mean)) # mean size of galls on a plant
diam.var = as.vector(by(data$gall_diameter, data$plant, var)) # mean size of galls on a plant
diam.median = as.vector(by(data$gall_diameter, data$plant, median)) # max size of galls on a plant
diam.max = as.vector(by(data$gall_diameter, data$plant, max)) # max size of galls on a plant

# tabulate gall diam just for healthy galls
gall.count.healthy = with(data[data$HEALTHY==1,], table(plant))
diam.mean.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, mean))) # mean size of galls on a plant
diam.median.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, median))) # median size of galls on a plant
diam.max.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, max))) # max size of galls on a plant

hist(diam.mean)
hist(diam.mean.healthy)

plot(galls ~ mean.per.4, data=data2)
ggplot(aes(x=mean.per.4, y= galls), data=data2)+
	geom_point()+
	stat_smooth()

#-------- plotting rates as a function of flws ---------------#

# survival rate as a func of # mean.per.4
ggplot(aes(x=mean.per.4, y=s), data=data2) +
	geom_point() +
	stat_smooth()

# prop of healthy to healthy+unk
ggplot(aes(x=mean.per.4, y=s2), data=data2) +
	geom_point() +
	stat_smooth()

# ptism rate as a func of # mean.per.4
ggplot(aes(x=mean.per.4, y=ptism.rate), data=data2) +
	geom_point() +
	stat_smooth()

# pred rate as a func of # mean.per.4
ggplot(aes(x=mean.per.4, y=pred.rate), data=data2) +
	geom_point() +
	stat_smooth()

# unk rate as a func of # mean.per.4
ggplot(aes(x=mean.per.4, y=unk.rate), data=data2) +
	geom_point() +
	stat_smooth()


## log scale #
# survival rate as a func of # log mean.per.4
ggplot(aes(x=log(mean.per.4+1), y=s), data=data2) +
	geom_point() +
	stat_smooth()
	
# prop of healthy to healthy+unk
ggplot(aes(x=log(mean.per.4+1), y=s2), data=data2) +
	geom_point() +
	stat_smooth()

# ptism rate as a func of # log mean.per.4
ggplot(aes(x=log(mean.per.4+1), y=ptism.rate), data=data2) +
	geom_point() +
	stat_smooth()

# pred rate as a func of # log mean.per.4
ggplot(aes(x=log(mean.per.4+1), y=pred.rate), data=data2) +
	geom_point() +
	stat_smooth()




#--------- gall size as a func of flws ---------------#

# mean diameter of healthy galls as a func of mean.per.4
ggplot(aes(x=mean.per.4, y=diam.mean.healthy),
	data=data2[data2$healthy>0,]) +
	geom_point() +
	stat_smooth()+
	stat_smooth(method='lm')

# max diameter of healthy galls as a func of mean.per.4
ggplot(aes(x=mean.per.4, y=diam.max),
	data=data2) +
	geom_point() +
	stat_smooth()+
	stat_smooth(method='lm')
	

#---------- models for mean gall size as func of flws -------#


m0 = lm(diam.mean.healthy[-28][-28] ~ 1)
m1 = lm(diam.mean.healthy[-28][-28] ~ data2[data2$healthy>0,'mean.per.4'][-28][-28])

AICctab(m0, m1, nobs=30)
AICtab(m0, m1)
anova(m0, m1)

plot(diam.mean.healthy[-28][-28] ~ data2[data2$healthy>0,][-28,]$mean.per.4[-28])
abline(m1)


a0 = lm(diam.mean[c(-22, -37)] ~ 1)
a1 = lm(diam.mean[c(-22,-37)] ~ data2$mean.per.4[c(-22,-37)])

AICctab(a0, a1, nobs=30)
anova(a0, a1)

plot(diam.mean[c(-22,-37)] ~ data2$mean.per.4[c(-22,-37)])
abline(a1)

