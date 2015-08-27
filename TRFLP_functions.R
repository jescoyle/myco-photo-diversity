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

	# If the curve is strictly positive, find the global minimum
	if(sum(use_cors<0)==0){
		local_mins = sapply(1:(length(divisors)-1), function(i){
			as.numeric(optimize(curve_func, c(divisors[i+1], divisors[i])))
		})
		use_div = local_mins[1,which(local_mins[2,]==min(local_mins[2,]))]
	}

	# If curve is strictly negative, find the global maximum
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





