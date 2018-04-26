# Fix fits

for (i in 1:81) {
	load(paste0('fit', i, '.RData'))
	fixedef <- subset(fit$coef, effect == 'fixed')
	randomef <- subset(fit$coef, effect == 'random')
	randomef <- randomef[,c(1,2,4,3,5)]
	names(randomef)[3:4] <- c('parameter','stat')
	fit$coef <- rbind(fixedef, randomef)
	save(fit, file = paste0('fit', i, '.RData'))
}

# Fix again

for (i in 1:81) {
	print(i)
	load(paste0('fit', i, '.RData'))
	fixedef <- fit$coef[1:36,1:5]
	randomef <- fit$coef[,6:10]
	fit$coef <- rbind(fixedef, randomef)
	save(fit, file = paste0('fit', i, '.RData'))
}
