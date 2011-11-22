##### Obs plots 2010-2011
# prelim analysis
# October 2011
library(ggplot2)
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

p1 = ggplot(data=d1, aes(x=galls2010, y=galls2011)) +
	geom_point() +
	geom_smooth(method='lm')
print(p1)

# let's stick with p2 for now since it has a larger sample size

p2 = ggplot(data=d2, aes(x=galls2010, y=galls2011)) +
	scale_x_continuous('Galls in 2010') + 
	scale_y_continuous('Galls in 2011') +
	geom_point(alpha = 0.35) +
	geom_smooth(method='lm', colour=1) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2)
ggsave(
	filename='~/Documents/Analysis repos/Obs plots/figs/plot2-time.ps',
	plot=p2, width=3.5, height=3)

m1 = lm(galls2011 ~ galls2010, data=d1)
m1 = mle2(galls2011 ~ dnorm(mean=a + b * galls2010, sd=sd), 
	start=list(a=2, b=1, sd=1), data=d1)
m0 = mle2(galls2011 ~ dnbinom(mu=a, size=s), 
	start=list(a=2, s=1), data=d1)
m2 = mle2(galls2011 ~ dnbinom(mu=a + b * galls2010, size=s), 
	start=list(a=2, b=1, s=1), data=d1)


AICtab(m2, m0)

# mean predictions and CI for m2
newg = 0:37
pred.m2 = data.frame(galls2010 = newg, 
	galls2011 = coef(m2)['a'] + coef(m2)['b'] * newg)

# confidence intervals for m2
# which values of r and a do I use for upper and lower bounds?
# I imagine, the bounds that give the lowest low and highest high?
profile.m2 = profile(m2)
ci.m2 = confint(profile.m2)
pred.m2$ymin = ci.m2['a', 1] + ci.m2['b', 1] * newg 
pred.m2$ymax = ci.m2['a', 2] + ci.m2['b', 2] * newg


p2 = ggplot(data=d2, aes(x=galls2010, y=galls2011)) +
	scale_x_continuous('Galls in 2010') + 
	scale_y_continuous('Galls in 2011') +
	geom_point(alpha = 0.35) +
	geom_smooth(aes(x=galls2010, y=galls2011, ymin=ymin, ymax=ymax), 
		data=pred.m2, stat='identity', colour=1, alpha=0.2) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2)
ggsave(
	filename='~/Documents/Analysis repos/Obs plots/figs/plot2-time.pdf',
	plot=p2, width=3.5, height=3)


p2space = ggplot(data=d2, aes(x=x, y=y, z=galls2011)) +
	geom_point(alpha=0.5) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	#stat_contour()
print(p2space)
#########################

# quality variables

p2size = ggplot(data=d2, aes(x=vol, y=galls2011)) +
	scale_x_continuous('plant foliage volume (cubic meters)') + 
	scale_y_continuous('mean number of galls') +
	geom_point(alpha=0.5)+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2size)

msize0 = mle2(galls2011 ~ dnbinom(mu = m, size=s),
	start=list(m=7, s=1), data=d2)
msize1 = mle2(galls2011 ~ dnbinom(mu = a + b * vol, size=s),
	start=list(a=7, b=0, s=1), data=d2)
msize2 = mle2(galls2011 ~ dnbinom(mu = a * vol * exp(-b * vol), size=s),
	start=list(a=1, b=1/30, s=1), data=d2)

# mean predictions and CI for msize1
newvol = seq(0,0.9,length=50)
pred.msize1 = data.frame(vol = newvol, 
	galls2011 = coef(msize1)['a'] + coef(msize1)['b'] * newvol)

# confidence intervals for m1
profile.msize1 = profile(msize1)
ci.msize1 = confint(profile.msize1)
pred.msize1$ymin = ci.msize1['a', 1] + ci.msize1['b', 1] * newvol 
pred.msize1$ymax = ci.msize1['a', 2] + ci.msize1['b', 2] * newvol


p2size = ggplot(data=d2, aes(x=vol, y=galls2011)) +
	scale_x_continuous('Plant foliage volume (cubic meters)') + 
	scale_y_continuous('Number of galls in 2011') +
	geom_point(alpha = 0.35) +
	geom_smooth(aes(x=vol, y=galls2011, ymin=ymin, ymax=ymax), 
		data=pred.msize1, stat='identity', colour=1, alpha=0.2) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2size)
ggsave(
	filename='~/Documents/Analysis repos/Obs plots/figs/plot2-size.pdf',
	plot=p2size, width=3.5, height=3)


####### size with gall density instead of abundance
AICtab(msize0, msize1, msize2)

msizec0 = mle2(density ~ dnorm(mean = m, sd=s),
	start=list(m=7, s=1), data=d2)
msizec1 = mle2(density ~ dnorm(mean = a + b * vol, sd=s),
	start=list(a=7, b=0, s=1), data=d2)
msizec2 = mle2(density ~ dnorm(mean = a * vol * exp(-b * vol), sd=s),
	start=list(a=1, b=1/30, s=1), data=d2)
msizec3 = mle2(density ~ dnorm(mean = a / (b + vol), sd=s),
	start=list(a=300, b=0.1, s=14), data=d2)
msizec4 = mle2(density ~ dnorm(mean = a * exp(-b * vol), sd=s),
	start=list(a=3000, b=2.5, s=14), data=d2, method='SANN')

summary(msizec4)

AICtab(msizec0, msizec1, msizec2, msizec3, msizec4)

# mean predictions and CI for msize1
newvol = seq(0,0.9,length=50)
pred.msizec3 = data.frame(vol = newvol, density = coef(msizec3)['a'] + 
	coef(msize1)['b'] * newvol)

# confidence intervals for m1
profile.msize1 = profile(msize1)
ci.msize1 = confint(profile.msize1)
pred.msize1$ymin = ci.msize1['a', 1] + ci.msize1['b', 1] * newvol 
pred.msize1$ymax = ci.msize1['a', 2] + ci.msize1['b', 2] * newvol


p2sizec = ggplot(data=d2, aes(x=vol, y=density)) +
	scale_x_continuous('Plant foliage volume (cubic meters)') + 
	scale_y_continuous('Number of galls per cubic meter') +
	geom_point(alpha=0.5)+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2sizec)
ggsave(
	filename='~/Documents/Analysis repos/Obs plots/figs/plot2-sizec.pdf',
	plot=p2sizec, width=3.5, height=3)