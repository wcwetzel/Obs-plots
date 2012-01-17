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

# mean.per.4 is total number of flowers 
# counted on all 4 stalks divided by 4
# so plants with fewer than 4 stalks
# have some "zero flower" flowering stalks
# mean.per.1 is mean flowers per stalks counted

par(mfrow=c(1,2))
plot(mgalls ~ mean.per.4, data=d1)
plot(mgalls ~ mean.per.1, data=d1)

par(mfrow=c(1,2))
plot(galls2011 ~ mean.per.4, data=d1)
plot(galls2011 ~ mean.per.1, data=d1)



