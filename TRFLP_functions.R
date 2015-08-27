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

## NEEDS FIX: (8/27/17) CAN'T HANDLE WHEN CORRELATION IS STRICTLY NEGATIVE- NEEDS TO MAXIMUM
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
			stddev = sum(these_peaks^2)/(length(these_peaks)-1)
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





