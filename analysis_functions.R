## This script contains functions for the analysis for the myco-photo diversity project


# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=NA)
	}

	if(!is.na(labels)){
		dZ = labels[length(labels)]-labels[1]
		dz = dY/dZ
		Yposition = barends[2] + dz*(labels-labels[1])

		text(barends[3]+dX*0.5, Yposition, round(labels,2), pos=4, cex=mycex)
		
		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		labels.round = round(labels, 2)
		
		## Determine how many characters away to place title
		digits = max(nchar(round(labels, 2))) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		no.decimal = sum(largest == floor(largest))>0 # Does one largest label lack a decimal?
			if(!no.decimal) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		no.negative = sum(largest >= 0)>0 # Does one largest label lack a negative sign?
			if(!no.negative) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
}

# A function that performs RDA on Hellinger-transformed community data matrix
#	comm = community data matrix
#	env = environmental data matrix with rows in same order as comm
#	binary = flag indicating where RDA should be performed on presence/absence data
calc_rda = function(comm, env, binary=F){

	# Convert to presence/absence if indicated
	if(binary) comm = comm > 0
	
	# Transform abundances
	comm_std = decostand(comm, 'hellinger')

	# Conduct constrained ordination
	ord = rda(comm_std, env)

	# Test for significance
	ord_test = anova(ord, nperm=9999)
	F = ord_test$F[1]
	P = ord_test$'Pr(>F)'[1]
	
	# Return the R-square
	R2 = RsquareAdj(ord)$r.squared
	R2adj = RsquareAdj(ord)$adj.r.squared

	data.frame(R2, R2adj, F, P)
}

# A function that calculates the correlation between compositional dissimilarity and environmental distance
# Standardization by total abundance is currently not implemented but could be easily added using vegan's decostand
#	comm = community data matrix
# 	metric = vector of dissimilarity methods allowed by vegdist function in vegan
#	env = environmental data matrix with rows in the same order as comm
#	binary = convert community matrix to presence/absence?	
calc_diss_cor = function(comm, metric, env, binary=F){
	if(binary) comm = comm > 0
	
	sapply(metric, function(m){
		comm_diss = vegdist(comm, method=m)
		env_dist = dist(env)
		cor(comm_diss, env_dist, use='pairwise.complete.obs')
	})
}


make_plot = function(xlim, ylim, xlab=NULL, ylab=NULL, cex=1){	
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)
	axis(1)
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	if(!is.null(xlab)) mtext(xlab, 1, 2.5, cex=cex)
	if(!is.null(ylab)) mtext(ylab, 2, 3, cex=cex)

}
