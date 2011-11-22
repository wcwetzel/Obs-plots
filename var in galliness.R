## variation in galliness
# 15 Nov 2011

library(ggplot2)
library(bbmle)
#d1 = read.csv('~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
d2 = read.csv('~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol

d2$varg = apply(data.frame(d2$galls2010, d2$galls2011), 1, var)
d2$sdg = apply(data.frame(d2$galls2010, d2$galls2011), 1, sd)


with(d2, cbind(galls2010, galls2011, sdg))
hist(d2$sdg)


plot(sdg ~ mgalls, data=d2)
plot(sdg ~ vol, data=d2)
abline(lm(sdg ~ vol, data=d2))
plot(sdg ~ fsh, data=d2)



