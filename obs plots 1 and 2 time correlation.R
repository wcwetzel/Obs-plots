####### plotting abundances through time in Plots 1 and 2 ###########
#### 5 May 2012 #####



library(ggplot2)
library(bbmle)
library(rethinking)


#------- plot 1 first ----------#

d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
d1 = d1[!is.na(d1$galls2010),]
d1$mgalls = with(d1, (galls2011 + galls2010)/2 )
d1$vol = with(d1, d1/100 * d2/100 * h/100) # cubic meters
d1$density = d1$galls2011/d1$vol


m0d1 = mle2(galls2011 ~ dnbinom(mu=a, size=s), 
	start=list(a=2, s=1), data=d1)
m1d1 = mle2(galls2011 ~ dnbinom(mu=a + b * galls2010, size=s), 
	start=list(a=2, b=1, s=1), data=d1)

AICtab(m1d1,m0d1)

post.m1d1 = sample.naive.posterior(m1d1)
new.2010 = seq(0,20, by=0.1)

m1d1.mu = sapply(new.2010,
	function(z) mean(post.m1d1[,1] + post.m1d1[,2] * z))

m1d1.ci = sapply(new.2010,
	function(z) HPDI(post.m1d1[,1] + post.m1d1[,2] * z))


p1NB = ggplot(data=d1, aes(x=galls2010, y=galls2011)) +
	scale_x_continuous('Galls in 2010') + 
	scale_y_continuous('Galls in 2011') +
	geom_point(alpha = 0.35) +
	geom_line(aes(x=new.2010, y=m1d1.mu)) +
	geom_line(aes(x=new.2010, y=m1d1.ci[1,]), lty=2) +
	geom_line(aes(x=new.2010, y=m1d1.ci[2,]), lty=2) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p1NB)
ggsave(
	filename='~/Documents/Analysis-repos/Obs-plots/figs/plot1NB-time.pdf',
	plot=p1NB, width=3.5, height=3)




#----------- second for plot 2 ----------------#


d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol


m0d2 = mle2(galls2011 ~ dnbinom(mu=a, size=s), 
	start=list(a=2, s=1), data=d2)
m1d2 = mle2(galls2011 ~ dnbinom(mu=a + b * galls2010, size=s), 
	start=list(a=2, b=1, s=1), data=d2)

AICtab(m1d2,m0d2)

post.m1d2 = sample.naive.posterior(m1d2)
new.2010 = seq(0,37, by=0.1)

m1d2.mu = sapply(new.2010,
	function(z) mean(post.m1d2[,1] + post.m1d2[,2] * z))

m1d2.ci = sapply(new.2010,
	function(z) HPDI(post.m1d2[,1] + post.m1d2[,2] * z))


p2NB = ggplot(data=d2, aes(x=galls2010, y=galls2011)) +
	scale_x_continuous('Galls in 2010') + 
	scale_y_continuous('Galls in 2011') +
	geom_point(alpha = 0.35) +
	geom_line(aes(x=new.2010, y=m1d2.mu)) +
	geom_line(aes(x=new.2010, y=m1d2.ci[1,]), lty=2) +
	geom_line(aes(x=new.2010, y=m1d2.ci[2,]), lty=2) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))
print(p2NB)
ggsave(
	filename='~/Documents/Analysis-repos/Obs-plots/figs/plot2NB-time.pdf',
	plot=p2NB, width=3.5, height=3)



