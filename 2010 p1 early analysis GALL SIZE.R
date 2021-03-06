# prelim look at 2010 p1 galls -- GALL SIZE
# 10 Apr 2012

library(ggplot2)
library(bbmle)
library(rethinking)

#------------------------- Data -----------------------#

### data by gall ###
data = read.csv('~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis.csv')
data = data[order(data$plant),] # sort data on plant ID number
# approx gall as a prolate spheroid:
data$gvol = (4/3) * pi * (data$gall_diameter/2)^2 * data$gall_length/2
# categorical variable for gall fate
data$fate[data$HEALTHY==1] = 'healthy'
data$fate[data$PARASITIZED==1] = 'parasitoid'
data$fate[data$PREDATION==1] = 'predator'
data$fate[data$UNK_MORT==1] = 'unk'



### data by gall, only for healthy or ptized galls ###
dph = data[data$HEALTHY == 1 | data$PARASITIZED == 1,]

### data by plant ###
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

#write.csv(data2, '~/Documents/DATA/2010 DATA/LAB/2010 plot1 gall dissections/2010 plot1 gall dissections WETZEL early analysis flws.csv')


#-------ploting histograms or densities of survival and attack rates by gall size by gall -------#


hist(data$gall_length)
hist(data$gall_diameter)

plot(HEALTHY ~ gvol, data=data)
plot(PARASITIZED ~ gvol, data=data)
plot(PREDATION ~ gvol, data=data)
plot(UNK_MORT ~ gvol, data=data)

par(mfrow=c(2,2))
plot(HEALTHY ~ gall_diameter, data=data)
plot(PARASITIZED ~ gall_diameter, data=data)
plot(PREDATION ~ gall_diameter, data=data)
plot(UNK_MORT ~ gall_diameter, data=data)

# this one is bad, cant compare between groups
p1 = ggplot(data=data, aes(x=gall_diameter, colour=fate, fill=fate, group=fate)) +
	geom_density(alpha=0.1, adjust=0.7, group=data$fate) + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p1

# this one is good!!
p1.count = ggplot(data=data, aes(x=gall_diameter, group=fate, colour = fate, fill=fate)) +
	scale_colour_brewer() +
	scale_fill_brewer() +
	geom_density(alpha=0.4, adjust=0.7, group=data$fate) + 
	aes(y=..count..) +
	scale_y_continuous('Count')+
	scale_x_continuous('Gall diameter (mm)') +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0, size=12),
	axis.title.y = theme_text(vjust = 0.5, angle=90, size=12))
p1.count

p1.fill = ggplot(data=data, aes(x=gall_diameter, colour=fate, fill=fate)) +
	geom_density(alpha=0.1, adjust=0.7, position='fill') + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p1.fill

# this one is good too!
p1.fillcount = ggplot(data=data, aes(x=gall_diameter, colour=fate, fill=fate)) +
	aes(y=..count..)+
	scale_colour_brewer() +
	scale_fill_brewer() +
	geom_density(adjust=0.7, position='fill') + theme_bw() +
	scale_y_continuous('Proportion')+
	scale_x_continuous('Gall diameter (mm)') +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0, size=12),
	axis.title.y = theme_text(vjust = 0.5, angle=90, size=12))
p1.fillcount


p1.fill.nopred = ggplot(data=data[data$fate != 'predator',], aes(x=gall_diameter, colour=fate, fill=fate)) +
	geom_density(alpha=0.1, adjust=0.7, position='fill') + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p1.fill.nopred

p2 = ggplot(data=data, aes(x=gall_diameter, colour=fate, fill=fate)) +
	geom_bar() + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p2

p3 = ggplot(data=data, aes(x=gall_diameter, colour=fate, fill=fate)) +
	geom_bar(position='fill') + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p3

p4 = ggplot(data=data[data$fate != 'predator',], aes(x=gall_diameter, colour=fate, fill=fate)) +
	geom_bar(position='fill') + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
p4


#--------------bin data and calculate precent parasitism, survival, predation for each bin------------#

data$diamBins = cut(data$gall_diameter, breaks=seq(3, 11, by=1))
with(data, table(diamBins, fate))
binData = reshape(as.data.frame(with(data, table(diamBins, fate))), v.names='Freq', idvar= 'diamBins',
	direction='wide', timevar='fate')
binData$diam = seq(3.5, 10.5, by=1)

# convert counts to proportions
binData$ptism = with(binData, Freq.parasitoid/(Freq.parasitoid + Freq.healthy))
binData$survivalFromPtism = with(binData, Freq.healthy/(Freq.parasitoid + Freq.healthy))

binData$survivalRate = with(binData, Freq.healthy/(Freq.parasitoid + Freq.healthy + Freq.predator +
	Freq.unk))
binData$predRate = with(binData, Freq.predator/(Freq.parasitoid + Freq.healthy + Freq.predator))
binData$ptismTot = with(binData, Freq.parasitoid/(Freq.parasitoid + Freq.healthy + Freq.predator +
	Freq.unk))
binData$predRateTot = with(binData, Freq.predator/(Freq.parasitoid + Freq.healthy + Freq.predator +
	Freq.unk))
binData$unkRate = with(binData, Freq.unk/(Freq.parasitoid + Freq.healthy + Freq.predator + Freq.unk))



# plot
plot(ptism ~ diamBins, data=binData, type='b')
plot(survivalRate ~ diamBins, data=binData)
plot(predRate ~ diamBins, data=binData)

# just ptism and healthy
plot(ptism ~ I(3:10), data=binData, type='b', col='red', ylab='rate', xlab='Binned gall diameter (mm)')
points(survivalFromPtism ~ I(3:10), data=binData, type='b', col='blue', pch=20)
legend(8, 1, c('Healthy', 'Parasitized'), col=c('blue', 'red'), lty=1, pch=c(21, 20))

# total
plot(ptismTot ~ seq(3.5,10.5,by=1), data=binData, type='b', col='red', ylab='rate', 
	xlab='Binned gall diameter (mm)', ylim=c(0,0.7))
points(survivalRate ~ seq(3.5,10.5,by=1), data=binData, type='b', col='blue', pch=20)
points(predRateTot ~ seq(3.5,10.5,by=1), data=binData, type='b', col='yellow', pch=23)
points(unkRate ~ seq(3.5,10.5,by=1), data=binData, type='b', col='violet', pch=15)
legend(8, 1, c('Healthy', 'Parasitized', 'Predated', 'Unknown'), col=c('blue', 'red', 'yellow', 'violet'), 
	lty=1, pch=c(20, 21, 23, 15))


ggplot(data=binData, aes(x=diamBins, y=ptism)) +
	geom_point() + stat_smooth(aes(group=1)) + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))

ggplot(data=binData, aes(x=diamBins, y=survivalRate)) +
	geom_point() + stat_smooth(aes(group=1)) + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))

ggplot(data=binData, aes(x=diamBins, y=predRate)) +
	geom_point() + stat_smooth(aes(group=1)) + theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))

#---------- diameter histogram plus binned survivorship and parasitism figure -------------#

par(mar=c(5,4,4,5)+0.1)
hist(data$gall_diameter, xlim=c(3,11), xlab='Gall diameter (mm)', 
	ylab='Frequency', main='')
par(new=TRUE)
plot(survivalRate ~ diam, data=binData, xlim=c(3,11), xaxt='n', yaxt='n', 
	xlab='', ylab='', type='b', col='green')
axis(4, las=1)
mtext('Rate', side=4, line=3)

par(new=TRUE)
plot(ptismTot ~ diam, data=binData, xlim=c(3,11), xaxt='n', yaxt='n', 
	xlab='', ylab='', type='b', col='blue')

par(new=TRUE)
plot(predRateTot ~ diam, data=binData, xlim=c(3,11), xaxt='n', yaxt='n', 
	xlab='', ylab='', type='b', col='red')

par(new=TRUE)
plot(unkRate ~ diam, data=binData, xlim=c(3,11), xaxt='n', yaxt='n', 
	xlab='', ylab='', type='b', col='orange')





#-------------- Binom modeling survival and attack probs by gall size ----------------------#
m0h = mle2(HEALTHY ~ dbinom(size=1, prob=plogis(x)), start=list(x=0), data=data)
m1h = mle2(HEALTHY ~ dbinom(size=1, prob=plogis(x + a * gall_length)), 
	start=list(x=0, a=0), data=data)
m2h = mle2(HEALTHY ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)), 
	start=list(x=0, a=0), data=data)
m3h = mle2(HEALTHY ~ dbinom(size=1, prob=plogis(x + a * gvol)), 
	start=list(x=0, a=0), data=data)

AICtab(m0h, m1h, m2h, m3h)


m0z = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x)), start=list(x=0), data=data)
m1z = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x + a * gall_length)), 
	start=list(x=0, a=0), data=data)
m2z = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)), 
	start=list(x=0, a=0), data=data)
m3z = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x + a * gvol)), 
	start=list(x=0, a=0), data=data)

AICtab(m0z, m1z, m2z, m3z)


m0p = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x)), start=list(x=0), data=data)
m1p = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x + a * gall_length)), 
	start=list(x=0, a=0), data=data)
m2p = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)), 
	start=list(x=0, a=0), data=data)
m3p = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x + a * gvol)), 
	start=list(x=0, a=0), data=data)

AICtab(m0p, m1p, m2p, m3p)


m0u = mle2(UNK_MORT ~ dbinom(size=1, prob=plogis(x)), start=list(x=0), data=data)
m1u = mle2(UNK_MORT ~ dbinom(size=1, prob=plogis(x + a * gall_length)), 
	start=list(x=0, a=0), data=data)
m2u = mle2(UNK_MORT ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)), 
	start=list(x=0, a=0), data=data)
m3u = mle2(UNK_MORT ~ dbinom(size=1, prob=plogis(x + a * gvol)), 
	start=list(x=0, a=0), data=data)

AICtab(m0u, m1u, m2u, m3u)

newdiam = seq(0, 12, length=1000)
pm2h = predict(m2h, newdata=data.frame(gall_diameter=newdiam))
pm2z = predict(m2z, newdata=data.frame(gall_diameter=newdiam))
pm0z = predict(m0z, newdata=data.frame(gall_diameter=newdiam))
pm2p = predict(m2p, newdata=data.frame(gall_diameter=newdiam))
pm2u = predict(m2u, newdata=data.frame(gall_diameter=newdiam))

par(mfrow=c(2,2))
plot(HEALTHY ~ gall_diameter, data=data)
lines(newdiam, pm2h)
plot(PARASITIZED ~ gall_diameter, data=data)
lines(newdiam, pm2z)
abline(h = pm0z)
plot(PREDATION ~ gall_diameter, data=data)
lines(newdiam, pm2p)
plot(UNK_MORT ~ gall_diameter, data=data)
lines(newdiam, pm2u)


#---------- maybe I need a multinomial model of some kind ------------#

x = with(data, cbind(HEALTHY, PARASITIZED, PREDATION, UNK_MORT))
with(data, cbind(x, apply(x, 1, sum), plant))
xt = t(x)

mnom = mle2(x ~ dmultinom(size=1, prob=c(h, z, p, u)), start=list(h=0.25,
	z=0.25, p=0.25, u=0.25), data=data)
	

#----------- how about just healthy vs. ptized? --------------------#

### data by gall, only for healthy or ptized galls ###
dph = data[data$HEALTHY == 1 | data$PARASITIZED == 1,]


plot(HEALTHY ~ gall_diameter, data=dph)
plot(PARASITIZED ~ gall_diameter, data=dph)

mph0 = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x)),
	start=list(x=0), data=dph)
mph1 = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)),
	start=list(x=0, a = 0), data=dph)
mph2 = mle2(PARASITIZED ~ dbinom(size=1, prob=plogis(x + a * gall_diameter + 
	b * gall_diameter^2)), start=list(x=0, a=0, b=0), data=dph)

AICtab(mph0, mph1, mph2)
anova(mph0, mph1)


# sample naive posterior for means and CIs
post.mph1 = sample.naive.posterior(mph1)
newdiam2 = seq(0, 12, length=1000)
mu.mph1 = sapply(newdiam2, function(z) mean(plogis(post.mph1[,1] + 
	post.mph1[,2] * z)))
ci.mph1 = sapply(newdiam2, function(z) HPDI(plogis(post.mph1[,1] + 
	post.mph1[,2] * z)))
	
# plot
means = with(dph, by(gall_diameter, PARASITIZED, mean))
medians = with(dph, by(gall_diameter, PARASITIZED, median))

par(mfrow=c(1,2))
plot(PARASITIZED ~ gall_diameter, data=dph, pch='|')
#points(medians, c(0,1), pch=18, col=3, cex=2)
lines(newdiam2, mu.mph1)
lines(newdiam2, ci.mph1[1,], lty=2)
lines(newdiam2, ci.mph1[2,], lty=2)


#----------- how about predation vs. everything except unkmort ---------------#
### data by gall, only for healthy or ptized galls ###
dnu = data[data$HEALTHY == 1 | data$PARASITIZED == 1 | data$PREDATION ==1,]

#plot(PREDATION ~ gall_diameter, data=dnu)


mpnu0 = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x)),
	start=list(x=0), data=dnu)
mpnu1 = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x + a * gall_diameter)),
	start=list(x=0, a = 0), data=dnu)
mpnu2 = mle2(PREDATION ~ dbinom(size=1, prob=plogis(x + a * gall_diameter + 
	b * gall_diameter^2)), start=list(x=0, a=0, b=0), data=dnu)

AICtab(mpnu0, mpnu1, mpnu2)
anova(mpnu0, mpnu1)


# sample naive posterior for means and CIs
post.mpnu1 = sample.naive.posterior(mpnu1)
newdiam3 = seq(0, 12, length=1000)
mu.mpnu1 = sapply(newdiam3, function(z) mean(plogis(post.mpnu1[,1] + 
	post.mpnu1[,2] * z)))
ci.mpnu1 = sapply(newdiam3, function(z) HPDI(plogis(post.mpnu1[,1] + 
	post.mpnu1[,2] * z)))
	
# plot
means2 = with(dnu, by(gall_diameter, PREDATION, mean))
medians2 = with(dnu, by(gall_diameter, PREDATION, median))

plot(PREDATION ~ gall_diameter, data=dnu, pch='|')
#points(medians2, c(0,1), pch=18, col=3, cex=2)
lines(newdiam3, mu.mpnu1)
lines(newdiam3, ci.mpnu1[1,], lty=2)
lines(newdiam3, ci.mpnu1[2,], lty=2)


#----- is var in gall size greater between plants than within plants? -------#

plot(gall_diameter ~ plant, data=data[data$plant<1300,])

data$fatecolor[data$fate=='healthy'] = 'green'
data$fatecolor[data$fate=='parasitoid'] = 'purple'
data$fatecolor[data$fate=='predator'] = 'orange'
data$fatecolor[data$fate=='unk'] = 'wheat'
plot(gall_diameter ~ plant, data=data[data$plant<1300,], pch=data$fate, col=data$fatecolor)

plot(gall_diameter ~ plant, data=data[data$plant<1300 & data$HEALTHY==1,])

plot(gall_diameter ~ as.factor(plant), data=data)
plot(gall_diameter ~ as.factor(plant), data=data[data$HEALTHY==1,])



a = aov(gall_diameter ~ as.factor(plant), data=data)
summary(a)

ah = aov(gall_diameter ~ as.factor(plant), data=data[data$HEALTHY==1,])
summary(ah)
