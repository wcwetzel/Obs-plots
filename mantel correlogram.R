### Mantel correlogram
# 6 Nov 2011

library(spatstat)
library(vegan)
library(compiler)


#d1 = read.csv(
#	'~/Documents/DATA/2010 DATA/FIELD/plot 1/plot1-2010-2011.csv')
#d2 = d1

d2 = read.csv(
	'~/Documents/DATA/2010 DATA/FIELD/plot 2/plot2-2010-2011.csv')
d2 = d2[!is.na(d2$galls2010),]
d2$mgalls = with(d2, (galls2011 + galls2010)/2 )
d2$vol = with(d2, d1/100 * d2/100 * h/100) # cubic meters
d2$density = d2$galls2011/d2$vol




n = nrow(d2)
dists = pairdist(d2[,c('x','y')])
p1 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=FALSE)
p2 = matrix(d2[,'galls2010'], nrow=n, ncol=n, byrow=TRUE)
diff = abs(p1 - p2)

comp.mc = cmpfun(mantel.correlog)
mc = comp.mc(D.eco = diff, D.geo = dists)

boot.mc = function(x){
	nboot = x
	outmat = matrix(nrow=16, ncol=nboot)
	
	for(i in 1:nboot){
		gboot = sample(d2$galls2010, replace=TRUE)
		p1 = matrix(gboot, nrow=n, ncol=n, byrow=FALSE)
		p2 = matrix(gboot, nrow=n, ncol=n, byrow=TRUE)
		dboot = abs(p1 - p2)
		mcboot = comp.mc(dboot, dists)
		outmat[,i] = mcboot$mantel.res[,'Mantel.cor']
	}
	return(outmat)
}

comp.boot.mc = cmpfun(boot.mc)

outmat = comp.boot.mc(500)

outmat.noNA = outmat[1:9,]

uq = function(x){quantile(x, 0.975)}
lq = function(x){quantile(x, 0.025)}
upper = apply(outmat.noNA,1,uq)
lower = apply(outmat.noNA,1,lq)

plot(mc$mantel.res[1:9, 'Mantel.cor'] ~ 
	mc$mantel.res[1:9, 'class.index'], 
	ylim=c(-0.03, 0.03), type='b',
	main='Plot 2 2010')
points(upper ~ mc$mantel.res[1:9,'class.index'], type='l')
points(lower ~ mc$mantel.res[1:9,'class.index'], type='l')

write.csv(cbind(mc$mantel.res[1:9,], outmat.noNA), file=
	'~/Documents/Analysis repos/Obs plots/bootMC-plot2-2010.csv')

