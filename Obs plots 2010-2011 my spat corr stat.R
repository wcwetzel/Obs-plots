library(ggplot2)
library(bbmle)

#d1 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
#d2 = d1


d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol


# let's stick with p2 for now since it has a larger sample size

# library(spatstat)

# # create window (landscape) in which to put points
# window = with( d2, owin(xrange=c(min(x), max(x)), yrange=c(min(y),
	# max(y))) )

# # create the ppp dataset with spatial coordinates 
# #	and plant data 'marks'
# markdat = with( d2, ppp(x, y, window=window, marks=data.frame(
	# galls2010=galls2010, galls2011=galls2011, 
	# vol=vol, density=density)) )

#x = dist(coords(markdat))


weirdDstat = function(x){ 
	# x must be a dataframe with 1st and 2nd column
	# holding x and y coordinates and 3 column
	# holding a mark (e.g., gall abundances)
	n = nrow(x)
	np = n * (n-1) / 2
	dists = pairdist(x[,1:2])
	p1 = matrix(x[,3], nrow=n, ncol=n, byrow=FALSE)
	p2 = matrix(x[,3], nrow=n, ncol=n, byrow=TRUE)
	sdev = (p1 - p2)^2
	diff = abs(p1 - p2)
	
	out = data.frame(D = as.vector(dists[lower.tri(dists)]),
		sdev = as.vector(sdev[lower.tri(sdev)]),
		diff = as.vector(diff[lower.tri(diff)]))
	
	r = 1/(np-1) * sum( (diff - mean(diff))/sd(diff) * 
		(dists - mean(dists))/sd(dists) )
	return(list(out, r))
}

outtrue = weirdDstat(data.frame(d2$x, d2$y, d2$galls2011))
cortrue = cor(outtrue$D, outtrue$diff)
plot(diff ~ D, data=outtrue)
ggplot(data=outtrue, aes(x=D, y=diff)) +
	geom_point() +
	stat_smooth(span=0.01)

# bootstrap it!

bs = 1000 # number of boot samples
corvec = vector(length=bs)

for(i in 1:bs){
	bootsample = sample(d2$galls2011, replace=TRUE)
	outtemp = weirdDstat(data.frame(d2$x, d2$y, bootsample))
	corvec[i] = cor(outtemp$D, outtemp$diff)
}

hist(corvec, breaks=20)
abline(v=cortrue)
quantile(corvec, c(0.025, 0.95, 0.975))
cortrue




####################
## Mantel from NCEAS pdf
########

mantel = function(x){ 
	# x must be a dataframe with 1st and 2nd column
	# holding x and y coordinates and 3 column
	# holding a mark (e.g., gall abundances)
	n = nrow(x)
	np = n * (n-1) / 2
	dists = pairdist(x[,1:2])
	p1 = matrix(x[,3], nrow=n, ncol=n, byrow=FALSE)
	p2 = matrix(x[,3], nrow=n, ncol=n, byrow=TRUE)
	diff = abs(p1 - p2)
	
	r = 1/(np-1) * sum( (diff - mean(diff))/sd(diff) * 
		(dists - mean(dists))/sd(dists) )
		
	return(r)
}

rtrue = mantel(data.frame(d2$x, d2$y, d2$galls2011))

# bootstrap it!

bs = 1000 # number of boot samples
rvec = vector(length=bs)

for(i in 1:bs){
	bootsample = sample(d2$galls2011, replace=TRUE)
	rvec[i] = mantel(data.frame(d2$x, d2$y, bootsample))
}

hist(rvec, breaks=20)
abline(v=rtrue)
quantile(rvec, c(0.025, 0.95, 0.975))
rtrue
fn = ecdf(rvec)
fn(rtrue)
