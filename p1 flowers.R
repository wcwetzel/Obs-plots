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

d1flws = d1[!is.na(d1$mean.per.4),]

# mean.per.4 is total number of flowers 
# counted on all 4 stalks divided by 4
# so plants with fewer than 4 stalks
# have some "zero flower" flowering stalks
# mean.per.1 is mean flowers per stalks counted

par(mfrow=c(1,2))
plot(mgalls ~ mean.per.4, data=d1flws)
plot(mgalls ~ mean.per.1, data= d1flws)

par(mfrow=c(1,2))
plot(galls2011 ~ mean.per.4, data= d1flws)
plot(galls2011 ~ mean.per.1, data= d1flws)

ggplot(data= d1flws, aes(x=mean.per.4, y=galls2011)) +
	geom_point()+
	stat_smooth(method='auto') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls')

m1 = mle2(galls2011 ~ dnbinom(mu=a + b * mean.per.4, size=k),
	start=list(a=5, b=0, k=1), data= d1flws)

m0 = mle2(galls2011 ~ dnbinom(mu=a, size=k),
	start=list(a=5, k=1), data= d1flws)

AICtab(m1, m0)

anova(m1, m0)

ggplot(data= d1flws, aes(x=mean.per.4, y=galls2011)) +
	geom_point()+
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	#stat_smooth(method='auto') +
	#stat_smooth(method='lm') +
	scale_x_continuous('Mean number of flowers per stalk') +
	scale_y_continuous('Number of galls') +
	geom_abline(intercept=coef(m1)['a'], slope=coef(m1)['b'])
	

