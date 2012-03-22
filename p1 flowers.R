### plot1 flowering stalks
### 17 Jan 2012


library(bbmle)
library(ggplot2)

d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')

d1 = d1[!is.na(d1$galls2010),]
d1$mgalls = with(d1, (galls2011 + galls2010)/2 )
d1$vol = with(d1, d1/100 * d2/100 * h/100) # cubic meters
d1$density2011 = d1$galls2011/d1$vol
d1$mdensity = d1$mgalls / d1$vol
d1$flwdens = d1$flws / d1$vol

d1flws = d1[!is.na(d1$mean.per.4),]

# mean.per.4 is total number of flowers 
# counted on all 4 stalks divided by 4
# so plants with fewer than 4 stalks
# have some "zero flower" flowering stalks
# mean.per.1 is mean flowers per stalks counted


## first plotting ##
par(mfrow=c(2,2))
plot(mgalls ~ mean.per.4, data=d1flws)
plot(mgalls ~ mean.per.1, data= d1flws)
plot(mgalls ~ inflono, data=d1flws)
plot(mgalls ~ flws, data= d1flws)

plot(log(mdensity) ~ log(flwdens), data=d1flws)

par(mfrow=c(2,2))
plot(galls2011 ~ mean.per.4, data= d1flws)
plot(galls2011 ~ mean.per.1, data= d1flws)
plot(galls2011 ~ inflono, data=d1flws)
plot(galls2011 ~ flws, data= d1flws)

ggplot(data= d1flws, aes(x=mean.per.4, y=galls2011)) +
	geom_point()+
	stat_smooth(method='auto') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls')

ggplot(data= d1flws, aes(x=inflono, y=galls2011)) +
	geom_point()+
	stat_smooth(method='auto') +
	scale_x_continuous('Number of infloresences') +
	scale_y_continuous('Number of galls')

ggplot(data= d1flws, aes(x=flws, y=galls2011)) +
	geom_point()+
	stat_smooth(method='auto') +
	scale_x_continuous('Estimated number of flowers') +
	scale_y_continuous('Number of galls')


d1flws = d1flws[d1flws$flws<60000,]
d1flws = d1flws[d1flws$flwdens < 780000,]

## models ##

# mean.per.4 and mean.per.1
m1 = mle2(galls2011 ~ dnbinom(mu=a + b * mean.per.4, size=k),
	start=list(a=5, b=0, k=1), data= d1flws)

m11 = mle2(galls2011 ~ dnbinom(mu=a + b * mean.per.1, size=k),
	start=list(a=5, b=0, k=1), data= d1flws)

m0 = mle2(galls2011 ~ dnbinom(mu=a, size=k),
	start=list(a=5, k=1), data= d1flws)

AICtab(m1, m0, m11)

anova(m1, m0)

ggplot(data= d1flws, aes(x=mean.per.4, y=galls2011)) +
	geom_point()+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	#stat_smooth(method='auto') +
	stat_smooth(method='lm') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls') 
	#geom_abline(intercept=coef(m1)['a'], slope=coef(m1)['b'])

ggplot(data= d1flws, aes(x=mean.per.1, y=galls2011)) +
	geom_point()+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	#stat_smooth(method='auto') +
	stat_smooth(method='lm') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls') +
	geom_abline(intercept=coef(m11)['a'], slope=coef(m11)['b'])
	


# total flws

m2 = mle2(galls2011 ~ dnbinom(mu=exp(a) + b * flws, size=k),
	start=list(a=1, b=0, k=1), data= d1flws, trace=TRUE, method='SANN')


AICtab(m2, m0, m1)

anova(m2, m0)

ggplot(data= d1flws, aes(x=flws, y=galls2011)) +
	geom_point()+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	stat_smooth(method='auto') +
	stat_smooth(method='lm') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls') +
	geom_abline(intercept=exp(coef(m2)['a']), slope=coef(m2)['b'])

# flws and volume

m3 = mle2(galls2011 ~ dnbinom(mu=exp(a) + b * flws + c * vol, size=k),
	start=list(a=1, b=0, c=0, k=1), data= d1flws, trace=TRUE, method='SANN')


AICtab(m3, m2, m1, m0)

anova(m3, m0)

plot(galls2011 ~ I(flws/vol), data=d1flws)

ggplot(data= d1flws, aes(x=flws, y=galls2011)) +
	geom_point()+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	#stat_smooth(method='auto') +
	#stat_smooth(method='lm') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls') +
	geom_abline(intercept=exp(coef(m3)['a']), slope=coef(m3)['b'])
