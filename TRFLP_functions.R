## This script contains functions for the analysis of TRFLP data


# A function that returns captured string from regexpr
extract_vals = function(text, regr){
	capstart = attr(regr, 'capture.start')
	capend = attr(regr, 'capture.length') + capstart - 1
	cap = data.frame(text = text, match = as.numeric(regr), start = as.numeric(capstart), stop = as.numeric(capend))
	
	ifelse(cap$match>0, substr(cap$text, cap$start, cap$stop), NA)	
}


# A function that thresholds peaks using constant percentage method of Sait et al 2003
# Chooses a baseline for all profiles that minimizes the correlation between number of peaks and total peak area
# peak_mat : a matrix of peak signal strength where rows are samples and columns are peak identities (e.g. size
thresh_pct_const = function(peak_mat){
	
	# Standardize each trace by its total signal strength
	tot_sig = apply(peak_mat, 1, sum)
	std_mat = peak_mat / tot_sig

	# Ordered vector of signal strengths
	signals = unique(std_mat[peak_mat>0])
	signals = signals[order(signals)]

	# For each signal strength, calculate the correlation between total peak area and num. peaks
	cors = sapply(signals, function(x){
		thresh_mat = peak_mat
		thresh_mat[std_mat<x] = 0
		npeaks = apply(thresh_mat, 1, function(y) sum(y>0))
		tot_sig = apply(thresh_mat, 1, sum)
		cor(tot_sig, npeaks)
	})

	# Which signal strength yields the smallest correlation?
	threshold = signals[which(abs(cors)==min(abs(cors), na.rm=T))][1]
	print(paste('Used threshold', round(threshold, 4), sep=' = '))

	# Make a plot
	plot(cors~signals, type='l', xlab='Proportion of total peak strength', ylab='Cor(Tot. peak strength, Num. peaks)')
	abline(v=threshold, col=2)
	abline(h=0)

	# Return thresholded data matrix
	peak_mat[std_mat<threshold] = 0
	peak_mat
	
}

# A function that thresholds peaks using the variable percentage threshold method of Osborne et al 2006
# Chooses a baseline for each profile that minimizes the correlation between number of peaks and total peak area
# peak_mat : a matrix of peak signal strength where rows are samples and columns are peak identities (e.g. size
thresh_pct_var = function(peak_mat){

	# Standardize each trace by its total signal strength
	tot_sig = apply(peak_mat, 1, sum)
	std_mat = peak_mat / tot_sig

	# Choose set of divisors
	div_start = signif(mean(tot_sig)*1000, 1)
	pow_start = floor(log10(div_start))
	pow_end = pow_start-2
	divisors = signif(10^seq(pow_start, pow_end, -.05), 2)

	# Calculate percentage thresholds in each profile based on divisors
	pcts = tot_sig %*% t(1/divisors)
	
	# Find the divisor that minimizes the correlation between number of peaks and total signal in each profile	
	cors = sapply(1:length(divisors), function(j){
		thresh_mat = peak_mat
		for(i in 1:nrow(thresh_mat)){
			thresh_mat[i, std_mat[i,] < pcts[i,j]] = 0
		}
		tot_area = rowSums(thresh_mat)
		npeaks = apply(thresh_mat, 1, function(y) sum(y>0))

		cor(tot_area, npeaks)
	})
	
	# If the first divisor leads to neg correlation, try a larger divisor
	while(cors[1] < 0){
		div_start = div_start*100
		pow_start = floor(log10(div_start))
		pow_end = pow_start-2
		divisors = signif(10^seq(pow_start, pow_end, -.05), 2)
		pcts = tot_sig %*% t(1/divisors)
		cors = sapply(1:length(divisors), function(j){
			thresh_mat = peak_mat
			for(i in 1:nrow(thresh_mat)){
				thresh_mat[i, std_mat[i,] < pcts[i,j]] = 0
			}
			tot_area = rowSums(thresh_mat)
			npeaks = apply(thresh_mat, 1, function(y) sum(y>0))

			cor(tot_area, npeaks)
		})
	}

	# Make a plot
	plot(cors~divisors, xlab='Divisor', ylab='Cor(Tot. peak strength, Num. peaks)')
	abline(h=0)

	# Which divisor yields the smallest correlation?
	curve_func = splinefun(divisors[!is.na(cors)], cors[!is.na(cors)])	
	curve(curve_func, add=T)

	# Ignore NA correlations
	use_cors = na.omit(cors)	

	# If the curve crosses 0, find the largest root
	if(sum(use_cors<0)>0&sum(use_cors>0)>0){
		lower_bound = max(divisors[cors<0], na.rm=T)
		upper_bound = divisors[1]
		use_div = uniroot(curve_func, c(lower_bound, upper_bound))$root
	
	} else {

	# Otherwise, If the curve is strictly positive, find the global minimum
	if(sum(use_cors<0)==0){
		local_mins = sapply(1:(length(divisors)-1), function(i){
			as.numeric(optimize(curve_func, c(divisors[i+1], divisors[i])))
		})
		use_div = local_mins[1,which(local_mins[2,]==min(local_mins[2,]))]
	}

	#Or, If curve is strictly negative, find the global maximum
	if(sum(use_cors>0)==0){
		local_maxs = sapply(1:(length(divisors)-1), function(i){
			as.numeric(optimize(curve_func, c(divisors[i+1], divisors[i]), maximum=T))
		})
		use_div = local_maxs[1,which(local_maxs[2,]==max(local_maxs[2,]))]
	}
	} # closes initial if-else statement


	abline(v=use_div, col=2)
	print(paste('Used divisor', use_div))

	# Calculate new percentages based on new divisor
	pcts = tot_sig %*% t(1/use_div)

	# Return thresholded data matrix
	for(i in 1:nrow(peak_mat)){
		peak_mat[i, std_mat[i,] < pcts[i,]] = 0
	}

	peak_mat
}


# A function that thresholds peaks using the statistical method of Adbo et al 2006.
# Peaks falling outside a predefined number of standard devitions are iteratively identified as true peaks.
# peak_mat : a matrix of peak signal strength where rows are samples and columns are peak identities (e.g. size
# sfact : the number of standard deviations used to identify peaks that differ from noise

thresh_stdev = function(peak_mat, sfact){
	
	# Standardize each trace by its total signal strength
	tot_sig = apply(peak_mat, 1, sum)
	std_mat = peak_mat / tot_sig

	# Iteratively identify peaks in the upper 3rd std dev of signal strength in the standardized data
	thresh_mat = std_mat
	for(i in 1:nrow(std_mat)){
		
		keep_peaks = rep(F, ncol(std_mat))
		names(keep_peaks) = colnames(std_mat)
		noise_peaks = std_mat[i,]
		foundpeaks = T	
		while(foundpeaks){

			these_peaks = noise_peaks[noise_peaks>0]
			stddev = sqrt(sum(these_peaks^2)/(length(these_peaks)-1))
			add_peaks = names(which(these_peaks>=sfact*stddev))

			if(length(add_peaks)>0){
				keep_peaks[add_peaks] = T
				noise_peaks[add_peaks] = 0
			} else { foundpeaks=F }
		}
		
		thresh_mat[i,!keep_peaks] = 0
	}

	thresh_mat
}

# A function that plots TRFLP profiles supplied as a matrix
# prof_mat = matrix of profiles
# blank = name(s) of samples(s) containing no DNA
# control = name of samples containinig control DNA
plot_profiles = function(prof_mat, blank=NA, control=NA){
	N = nrow(prof_mat)

	par(mfrow=c(N,1))
	par(mar=c(2,4,0.5,1))
	if(!is.na(blank[1])){
		
		this_pf = prof_mat[blank,]
		if(length(blank)>1){
			this_pf = apply(this_pf, 2, function(x) x[x==max(x)][1])		
		}
		
		plot(as.numeric(names(this_pf)), this_pf, type='h', 
			ylim=c(0, max(this_pf)+1), ylab='Blank', xlab='', axes=F, col='blue')
		axis(2, las=1)
		blank_lines = names(this_pf[this_pf>0])
		prof_mat = prof_mat[!(rownames(prof_mat) %in% blank),]
	}
	if(!is.na(control)){
		this_pf = prof_mat[control,]
		plot(as.numeric(names(this_pf)), this_pf, type='h', 
			ylim=c(0, max(this_pf)+1), ylab='Control', xlab='', las=1, axes=F, col='red')
		axis(2, las=1)
		if(!is.na(blank)) points(blank_lines, rep(max(this_pf), length(blank_lines)), pch=3, col='blue')

		control_lines = as.numeric(names(this_pf[this_pf>0]))
		prof_mat = prof_mat[rownames(prof_mat)!=control,]
	}

	for(i in 1:nrow(prof_mat)){
		this_pf = prof_mat[i,]
		plot(as.numeric(names(this_pf)), this_pf, type='h', 
			ylim=c(0, max(this_pf)+1), ylab=rownames(prof_mat)[i], xlab='', las=1, axes=F)
		axis(2, las=1)
		if(i==nrow(prof_mat)) axis(1)
		
		if(!is.na(blank)) points(blank_lines, rep(max(this_pf), length(blank_lines)), pch=3, col='blue')
		if(!is.na(control)) points(control_lines, rep(max(this_pf), length(control_lines)), pch=3, col='red')
	
	}
}


## This function bins peaks that likely come from the same fragment length using hierarchical clustering

align_peaks = function(prof_mat, drawplot=T){
	#require(cluster)
	sizes = as.numeric(colnames(prof_mat))

	# Calculate distance matrix
	dmat = dist(sizes)

	# UPGMA clustering
	cl = hclust(dmat, method='average')
	
	# Define groups based on a distance of 1bp
	# NEED TO CHECK PUB FOR WHETHER THIS IS THE CORRECT WAY TO DO THIS
	# MAY WANT TO USE A LARGER WINDOW
	bins = cutree(cl, h=1)

	# Sum peak heights within bins
	new_mat = matrix(0, nrow=nrow(prof_mat), ncol=max(bins))
	for(i in 1:ncol(new_mat)){
		merge_cols = as.character(sizes[bins==i])
		if(length(merge_cols)>1){
			new_mat[,i] = rowSums(prof_mat[,merge_cols])
		} else {
			new_mat[,i] = prof_mat[,merge_cols]	
		}
	}

	# Determine fragment length for each bin
	new_sizes = sapply(1:max(bins), function(i) round(mean(sizes[bins==i]),0))
	colnames(new_mat) = new_sizes
	rownames(new_mat) = rownames(prof_mat)

	# Plot the clustering
	if(drawplot) plot(cl, labels=new_sizes[bins])

	# Remove columns that no longer have peaks
	new_mat =  new_mat[,colSums(new_mat)!=0]

	new_mat
}


# A function that drops peaks associated with a control
# control = name of row in prof_mat corresponding to control
drop_peaks = function(prof_mat, control){
	control_peaks = prof_mat[control,]>0
	prof_mat = as.matrix(prof_mat[,!control_peaks])
	prof_mat = as.matrix(prof_mat[rownames(prof_mat)!=control,])

	print(paste('Dropped', sum(control_peaks),'/',length(control_peaks), 'peaks'))
	prof_mat
}

# A function that tests the extent to which a given profile exists in another profile and th
# Names of the two profiles should be in bp
# pattern = profile to search for
# target = profile to search in
match_peaks = function(pattern, target){
	pattern = pattern[,pattern!=0]
	pattern_rank = rank(pattern)
	found_peaks = which(target[,names(pattern)]>0)
	if(length(found_peaks)>0){
		matched = names(pattern)[found_peaks]
		Pct_match = sum(pattern[,matched])/sum(pattern)
		N_match = length(matched)
	} else {
		Pct_match = 0
		N_match = 0
	}
	N_peaks = length(pattern)

	c(Pct_match=Pct_match, N_match=N_match, N_peaks=N_peaks)
}


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





