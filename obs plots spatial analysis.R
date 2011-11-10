##### Obs plots 2010-2011
# marked point process analysis
#### !!!!! ACTUALLY I DONT THINK THIS SHOULD BE TREATED AS A POINT
#### PROCESS B/C I WANT THE POSITIONS OF PLANTS TO BE FIXED NOT
#### RANDOM
# November 2011
library(ggplot2)
library(bbmle)

d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
d2 = d1


#d2 = read.csv(
#	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol


# let's stick with p2 for now since it has a larger sample size

library(spatstat)

# create window (landscape) in which to put points
window = with( d2, owin(xrange=c(min(x), max(x)), yrange=c(min(y),
	max(y))) )
# ppp dataset without marks
dat = with(d2, ppp(x, y, window=window))
# create the ppp dataset with spatial coordinates 
#	and plant data 'marks'
markdat = with( d2, ppp(x, y, window=window, marks=data.frame(
	galls2010=galls2010, galls2011=galls2011, vol=vol, density=density)) )

# 1st analyze spatial point pattern of Artr
# Ripley's K
ripK = Kest(dat)
plot(ripK)
ripK.ci = envelope(dat, 'Kest')
plot(ripK.ci)
# Diggle's G function (empty space)
G = Gest(dat)
G.ci = envelope(dat, 'Gest')
plot(G)
plot(G.ci)
# Diggle's F function (nearest neighbor)
F = Fest(dat)
F.ci = envelope(dat, 'Fest')
plot(F)
plot(F.ci)
# J function
J = Jest(dat)
J.ci = envelope(dat, 'Jest')
plot(J)
plot(J.ci)


# mark correlation function- looks at correlation among marks on events
# each spatial location is an event (a plant) and each mark is something
# about that event
mc1 = markcorr(markdat)
#plot(mc1)
#plot(mc1$galls2010$trans ~ mc1$galls2010$r, type='l')
mc2 = markcorr(markdat, f = function(m1, m2){(m1 - m2)^2})
#plot(mc2)


# BOOTSTRAP NULL MODEL
# reshuffle number of galls per plant
# but keep plant locations identical

nboot = 5000 # number of bootstrap replicates
boottrans2010 = matrix(ncol=nboot, nrow=length(mc1$galls2010$trans))
bootiso2010 = matrix(ncol=nboot, nrow=length(mc1$galls2010$trans))
boottrans2011 = matrix(ncol=nboot, nrow=length(mc1$galls2010$trans))
bootiso2011 = matrix(ncol=nboot, nrow=length(mc1$galls2010$trans))
bootdata = with(d2, ppp(x, y, window=window))

for(i in 1:nboot){
	bootgalls2010 = sample(d2$galls2010, size=nrow(d2), replace=TRUE)
	bootgalls2011 = sample(d2$galls2011, size=nrow(d2), replace=TRUE)
	marks(bootdata) = data.frame(bootgalls2010, bootgalls2011)
	mctemp = markcorr(bootdata)
	boottrans2010[,i] = mctemp$bootgalls2010$trans
	bootiso2010[,i] = mctemp$bootgalls2010$iso
	boottrans2011[,i] = mctemp$bootgalls2011$trans
	bootiso2011[,i] = mctemp$bootgalls2011$iso
}

# quantile functions with 0.025 and 0.975 as default
# for use with apply
quantile025 = function(x){quantile(x, probs=0.025)}
quantile975 = function(x){quantile(x, probs=0.975)}

# calculate upper and lower CIs for bootstrapped data
# at each distance value
uppertrans2010 = apply(boottrans2010, MARGIN=1, FUN=quantile975)
lowertrans2010 = apply(boottrans2010, MARGIN=1, FUN=quantile025)
upperiso2010 = apply(bootiso2010, MARGIN=1, FUN=quantile975)
loweriso2010 = apply(bootiso2010, MARGIN=1, FUN=quantile025)

uppertrans2011 = apply(boottrans2011, MARGIN=1, FUN=quantile975)
lowertrans2011 = apply(boottrans2011, MARGIN=1, FUN=quantile025)
upperiso2011 = apply(bootiso2011, MARGIN=1, FUN=quantile975)
loweriso2011 = apply(bootiso2011, MARGIN=1, FUN=quantile025)


# plot trans corrected 2010
par(mfrow=c(2,2))
plot(mc1$galls2010$trans ~ mctemp$bootgalls2010$r, type='l', 
	ylim=c(min(lowertrans2010), max(uppertrans2010)), las=1,
	main='PLOT 1: Galls2010 trans', xlab='Distance (r)', ylab='K(r)',
	sub='5000 bootstrap replicates')
points(uppertrans2010 ~ mctemp$bootgalls2010$r, type='l', lty=2)
points(lowertrans2010 ~ mctemp$bootgalls2010$r, type='l', lty=2)
abline(h=1, lty=3)

# plot iso corrected 2010
plot(mc1$galls2010$iso ~ mctemp$bootgalls2010$r, type='l', 
	ylim=c(min(loweriso2010), max(upperiso2010)), las=1,
	main='Galls2010 iso', xlab='Distance (r)', ylab='K(r)')
points(upperiso2010 ~ mctemp$bootgalls2010$r, type='l', lty=2)
points(loweriso2010 ~ mctemp$bootgalls2010$r, type='l', lty=2)
abline(h=1, lty=3)

# plot trans corrected 2011
plot(mc1$galls2011$trans ~ mctemp$bootgalls2011$r, type='l', 
	ylim=c(min(lowertrans2011), max(uppertrans2011)), las=1,
	main='Galls2011 trans', xlab='Distance (r)', ylab='K(r)')
points(uppertrans2011 ~ mctemp$bootgalls2011$r, type='l', lty=2)
points(lowertrans2011 ~ mctemp$bootgalls2011$r, type='l', lty=2)
abline(h=1, lty=3)

# plot iso corrected 2011
plot(mc1$galls2011$iso ~ mctemp$bootgalls2011$r, type='l', 
	ylim=c(min(loweriso2011), max(upperiso2011)), las=1,
	main='Galls2011 iso', xlab='Distance (r)', ylab='K(r)')
points(upperiso2011 ~ mctemp$bootgalls2011$r, type='l', lty=2)
points(loweriso2011 ~ mctemp$bootgalls2011$r, type='l', lty=2)
abline(h=1, lty=3)



