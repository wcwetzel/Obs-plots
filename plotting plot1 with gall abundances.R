logistic = function(x){1 / (1 + exp(x))}

plot(d1$y, d1$x, pch=20, ylab='meters', xlab='meters', 
	col=grey(trans(sqrt(d1$galls2011))))
points(d1$y, d1$x)
	
	

plot(trans(sqrt(galls2011)) ~ galls2011, data=d1)

trans = function(x) { (max(x) - x) / max(x)}

trans(d1$galls2011)
trans(sqrt(d1$galls2011))


max(d1$galls2011)