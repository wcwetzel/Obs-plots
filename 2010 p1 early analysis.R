# prelim look at 2010 p1 galls
# 7 Apr 2012


### some early results ######
# survival rate doesn't change or has a very slight decrease
# ptism rate doesn't vary with gall abundance, or maybe has a weak
# predation rate seems to increase slightly with gall abundance

# mean gall diameter doesn't increase with galliness
# mean diam of healthy galls doesn't increase with galliness
# max gall diameter does seem to have a strong increasing relationship 
	# with galliness, but I think this is just a sampling effect.
######

library(ggplot2)
library(bbmle)
library(rethinking)
library(lme4)


data = read.csv('~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis.csv')
data = data[order(data$plant),] # sort data on plant ID number

# tabulate gall outcomes by plant
gall.count = table(data$plant)
healthy.count = as.vector(by(data$HEALTHY, data$plant, sum))
ptoid.count = as.vector(by(data$PARASITIZED, data$plant, sum))
pred.count = as.vector(by(data$PREDATION, data$plant, sum))
unk.count = as.vector(by(data$UNK_MORT, data$plant, sum))
diam.mean = as.vector(by(data$gall_diameter, data$plant, mean)) # mean size of galls on a plant
diam.var = as.vector(by(data$gall_diameter, data$plant, var)) # mean size of galls on a plant
diam.median = as.vector(by(data$gall_diameter, data$plant, median)) # max size of galls on a plant
diam.max = as.vector(by(data$gall_diameter, data$plant, max)) # max size of galls on a plant


# tabulate gall diam just for healthy galls
gall.count.healthy = with(data[data$HEALTHY==1,], table(plant))
diam.mean.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, mean))) # mean size of galls on a plant
diam.median.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, median))) # median size of galls on a plant
diam.max.healthy = with(data[data$HEALTHY==1,], as.vector(by(gall_diameter, plant, max))) # max size of galls on a plant


# make new data frame in which each row is a plant
data2 = data.frame(
	plant = dimnames(gall.count)[[1]], 
	galls = as.vector(gall.count),
	healthy = healthy.count,
	ptism = ptoid.count,
	pred = pred.count,
	unk = unk.count,
	diam.mean=diam.mean,
	diam.median = diam.median,
	diam.max=diam.max
	)

write.csv(data2, "/Users/will/Desktop/plot1 with fates.csv")

# data frame for healthy galls
data2healthy = data.frame(
	plant = dimnames(gall.count.healthy)[[1]],
	galls = as.vector(gall.count)[data2$healthy>0],
	galls.healthy = as.vector(gall.count.healthy),
	diam.mean.healthy = diam.mean.healthy,
	diam.median.healthy = diam.median.healthy,
	diam.max.healthy = diam.max.healthy
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


#--------------------- plotting rates ----------------------------#

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

# unk rate as a func of # galls
ggplot(aes(x=galls, y=unk.rate), data=data2) +
	geom_point() +
	stat_smooth()


## log scale #
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





panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use='pairwise.complete.obs'))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 1) #cex.cor * r *3)
}

pairs(~ s + ptism.rate + pred.rate + unk.rate + galls + log(galls),
	data=data2, upper.panel=panel.smooth, lower.panel=panel.cor)






#----- is there a relationship between ptism and unkmort? ---------#

# by rate
ggplot(aes(x=ptism.rate, y=unk.rate), data=data2) +
	geom_point(position=position_jitter(w = 0.02, h = 0.02), alpha=.7) +
	scale_x_continuous('Parasitism rate') +
	scale_y_continuous('Unk mortality rate') +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0, size=12))


# by number of galls
ggplot(aes(x=ptism, y=unk), data=data2) +
	geom_point(position=position_jitter(w = 0.1, h = 0.1), alpha=.7) +
	scale_x_continuous('Parasitized galls') +
	scale_y_continuous('Unk mortality galls') +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0, size=12))

with(data2, cor(ptism, unk))
with(data2, cor(ptism.rate, unk.rate))


#--------------- some models of fates ~ galliness ---------------#


# !!!!!!! Should there be a mixed effect for plant? !!!!!!! #
# !!!!!!!!!!!!!!!!! YES !!!!!!!!!!!!!!!!!!! #

## survival ##
m0.healthy = mle2(healthy ~ dbinom(size=galls, prob=plogis(p)), data=data2, 
	start=list(p=0))
m0.healthy = lmer(cbind())

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
#points(predicted.p0 ~ galls, data=data2, pch='0', col='red')
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

AICctab(m0.pred, m1.pred, m2.pred, m3.pred, weights=TRUE, nobs=40)
anova(m0.pred, m1.pred)

newdata = data.frame(galls=seq(1,19,length=100))

predicted.p0 = predict(m0.pred,newdata=newdata)/newdata$galls
predicted.p1 = predict(m1.pred, newdata)/newdata$galls
predicted.p2 = predict(m2.pred, newdata)/newdata$galls
predicted.p3 = predict(m3.pred, newdata)/newdata$galls

post.m1.pred = sample.naive.posterior(m1.pred)

m1.pred.ci = sapply(newdata$galls,
	function(z) HPDI(plogis(post.m1.pred[,1] + post.m1.pred[,2] * z)))

plot(pred.rate + runif(40, -0.03,0.03) ~ I(galls + runif(40, -0.2,0.2)), data=data2, pch=20, las=1,
	ylab='Predation rate', xlab='Number of galls', type='n')
#points(predicted.p0 ~ galls, data=newdata, type='l', col='red')
#points(predicted.p1 ~ galls, data=newdata, type='l')
#points(predicted.p2 ~ galls, data=newdata, type='l', col='green')
#points(predicted.p3 ~ galls, data=newdata, type='l', col='yellow')
#lines(x=newdata$galls, y=m1.pred.ci[1,], lty=2)
#lines(x=newdata$galls, y=m1.pred.ci[2,], lty=2)
polygon(c(newdata$galls, rev(newdata$galls)), 
	c(m1.pred.ci[1,], rev(m1.pred.ci[2,])), col='grey', border=FALSE)
points(predicted.p1 ~ galls, data=newdata, type='l')
points(pred.rate + runif(40, -0.03,0.03) ~ I(galls + runif(40, -0.2,0.2)), data=data2, pch=20)

#------- plot a grid of all relationships -------------------#

par(mfcol = c(2, 4))

plot(s + runif(40, -0.03,0.03) ~ I(galls + runif(40, -0.2,0.2)), data=data2, 
	pch=20, ylab='Survival rate', xlab='Galls')
lines(lowess(data2$s ~ data2$galls), type='l')
plot(s + runif(40, -0.03,0.03) ~ log(galls + runif(40, 0, 0.05)), data=data2, 
	pch=20,  ylab='Survival rate', xlab='log(Galls)')
lines(lowess(data2$s ~ log(data2$galls)), type='l')

plot(ptism.rate + runif(40, -0.03,0.03)  ~ I(galls + runif(40, -0.2,0.2)), 
	data=data2, pch=20, ylab='Parasitism rate', xlab='Galls')
lines(lowess(data2$ptism.rate ~ data2$galls), type='l')
plot(ptism.rate + runif(40, -0.03,0.03) ~ log(galls + runif(40, 0, 0.05)), 
	data=data2, pch=20, ylab='Parasitism rate', xlab='log(Galls)')
lines(lowess(data2$ptism.rate ~ log(data2$galls)), type='l')

plot(pred.rate + runif(40, -0.03,0.03) ~ I(galls + runif(40, -0.2,0.2)), data=data2, pch=20, las=1,
	ylab='Predation rate', xlab='Galls')
lines(lowess(data2$pred.rate ~ data2$galls), type='l')
plot(pred.rate + runif(40, -0.03,0.03) ~ log(galls + runif(40, 0, 0.05)), 
	data=data2, pch=20,  ylab='Predation rate', xlab='log(Galls)')
lines(lowess(data2$pred.rate ~ log(data2$galls)), type='l')

plot(unk.rate + runif(40, -0.03,0.03) ~ I(galls + runif(40, -0.2,0.2)), data=data2, 
	pch=20, ylab='Unk mortality rate', xlab='Galls')
lines(lowess(data2$unk.rate ~ data2$galls), type='l')
plot(unk.rate + runif(40, -0.03,0.03) ~ log(galls + runif(40, 0, 0.05)), 
	data=data2, pch=20,  ylab='Unk mortality rate', xlab='log(Galls)')
lines(lowess(data2$unk.rate ~ log(data2$galls)), type='l')




#--------------- plotting gall diameter as a func of galliness --------------------#
# mean gall diameter as func of # galls
ggplot(aes(x=galls, y=diam.mean), data=data2) +
	geom_point() +
	stat_smooth()

# mean gall diameter of healthy galls as func of # galls
ggplot(aes(x=galls, y=diam.mean.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# mean gall diameter of healthy galls as func of # healthy galls
ggplot(aes(x=galls.healthy, y=diam.mean.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# median gall diameter as func of # galls
ggplot(aes(x=galls, y=diam.median), data=data2) +
	geom_point() +
	stat_smooth()

# median gall diameter of healthy galls as func of # galls
ggplot(aes(x=galls, y=diam.median.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# median gall diameter of healthy galls as func of # healthy galls
ggplot(aes(x=galls.healthy, y=diam.median.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# max gall diameter as func of # galls
ggplot(aes(x=galls, y=diam.max), data=data2) +
	geom_point() +
	stat_smooth()

# max gall diameter of healthy galls as func of # galls
ggplot(aes(x=galls, y=diam.max.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

## log scale ##
# mean gall diameter as func of # galls
ggplot(aes(x=log(galls), y=diam.mean), data=data2) +
	geom_point() +
	stat_smooth()

# mean gall diameter of healthy galls as func of # galls
ggplot(aes(x=log(galls), y=diam.mean.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# mean gall diameter of healthy galls as func of # healthy galls
ggplot(aes(x=log(galls.healthy), y=diam.mean.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# median gall diameter as func of # galls
ggplot(aes(x=log(galls), y=diam.median), data=data2) +
	geom_point() +
	stat_smooth()

# median gall diameter of healthy galls as func of # galls
ggplot(aes(x=log(galls), y=diam.median.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# median gall diameter of healthy galls as func of # healthy galls
ggplot(aes(x=log(galls.healthy), y=diam.median.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

# max gall diameter as func of # galls
ggplot(aes(x=log(galls), y=diam.max), data=data2) +
	geom_point() +
	stat_smooth()

# max gall diameter of healthy galls as func of # galls
ggplot(aes(x=log(galls), y=diam.max.healthy), data=data2healthy) +
	geom_point() +
	stat_smooth()

#------------------ modeling gall size as a function of galliness --------------#

s0 = lm(diam.mean.healthy ~ 1, data=data2healthy)
s1 = lm(diam.mean.healthy ~ galls, data=data2healthy)
s2 = lm(diam.mean.healthy ~ log(galls), data=data2healthy)

AICtab(s0, s1, s2)

sm0 = lm(diam.max ~ 1, data=data2)
sm1 = lm(diam.max ~ galls, data=data2)
sm2 = lm(diam.max ~ log(galls), data=data2)

AICtab(sm0, sm1, sm2)

post.sm2 = sample.naive.posterior(sm2)

newloggalls = seq(0, 3, length=1000)

sm2.mu = sapply(newloggalls,
	function(z) mean(post.sm2[,1] + post.sm2[,2] * z))

sm2.ci = sapply(newloggalls,
	function(z) HPDI(post.sm2[,1] + post.sm2[,2] * z))


plot(diam.max ~ log(galls), data=data2, pch=20, ylab="Largest diameter")
lines(newloggalls, sm2.mu)
lines(newloggalls, sm2.ci[1,], lty=2)
lines(newloggalls, sm2.ci[2,], lty=2)

