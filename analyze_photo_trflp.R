## This script is used to analyze photobiont T-RFLP profiles generated from the Photobiont-Mycobiont Diversity Project

options(stringsAsFactors=F)

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/TRFLP/photo_TRFLP/'
setwd(working_dir)

# Read in functions for the analysis of TRFLP data
source('TRFLP_functions.R')


# Read in peakscanner file and re-name columns
profiles_file = 'photo_trflp.txt'
peaksc_raw = read.table(paste(data_dir,profiles_file, sep=''), sep='\t', header=T, check.names=F)
peaks = peaksc_raw[,c('Sample Name','Dye/Sample Peak','Size','Height','Width in BP', 'Area in BP','UD3')]
names(peaks) = c('Sample','DyeNum','Size','Height','Width', 'Area','strain')

# Split Dye and Peak numbers
peaks$Dye = substr(peaks$DyeNum, 1, 1)
peaks$Num = as.numeric(extract_vals(peaks$DyeNum, regexpr('[A-Z], ([0-9]+)', peaks$DyeNum, perl=T)))


### Standardization, thresholding, and alignment

# Remove missing peaks
peaks = subset(peaks, !is.na(Size))

# Remove peaks less than 25 bp, which could be formed by primer dimers or remaining fluorophores
peaks = subset(peaks, Size > 25)

# Convert peak sizes to nearest integer
peaks$Size = round(peaks$Size, 0)

## Convert dataframe to array of abundance matrices based on peak area and peak height

# Define dyes to focus on
good_dye = c('B','Y')

# Tabulate peaks
height_mat = xtabs(Height~Sample+Size+Dye, data=subset(peaks, Dye %in% good_dye))
area_mat = xtabs(Area~Sample+Size+Dye, data=subset(peaks, Dye %in% good_dye))

# Examine correlation between total signal and number of peaks in each sample
# This determines whether thresholding needs to account for variations in amount of DNA
use_mat = area_mat[,,'B']
use_mat = use_mat[rownames(use_mat)!='empty',]
npeaks = apply(use_mat, 1, function(x) sum(x>0))
totsignal = apply(use_mat, 1, sum)
plot(totsignal, npeaks)
cor(totsignal, npeaks)

# Compare constant and variable percentage thresholding methods
const_pct_mat = thresh_pct_const(use_mat)
var_pct_mat = thresh_pct_var(use_mat)
par(mfrow=c(1,3))
npeaks = apply(use_mat, 1, function(x) sum(x>0))
totsignal = apply(use_mat, 1, sum)
plot(totsignal, npeaks)
npeaks = apply(const_pct_mat, 1, function(x) sum(x>0))
totsignal = apply(const_pct_mat, 1, sum)
plot(totsignal, npeaks)
npeaks = apply(var_pct_mat, 1, function(x) sum(x>0))
totsignal = apply(var_pct_mat, 1, sum)
plot(totsignal, npeaks)

# Statistical thresholding method does not work well with this data.
par(mfrow=c(2,5))
for(i in seq(3, 21,2)){
	thresh_mat = thresh_stdev(use_mat, i) 
	npeaks = apply(thresh_mat, 1, function(x) sum(x>0))
	totsignal = apply(var_pct_mat, 1, sum)
	plot(totsignal, npeaks, main=paste('Sfact = ', i))
	print(cor(totsignal, npeaks))
	abline(lm(npeaks~totsignal))
}


## Standardize abundance matrices by total area or total height
height_tot = apply(height_mat, c(1,3), sum)
height_mat[,,good_dye[1]] = height_mat[,,good_dye[1]]/height_tot[,good_dye[1]]
height_mat[,,good_dye[2]] = height_mat[,,good_dye[2]]/height_tot[,good_dye[2]]

area_tot = apply(area_mat, c(1,3), sum)
area_mat[,,good_dye[1]] = area_mat[,,good_dye[1]]/area_tot[,good_dye[1]]
area_mat[,,good_dye[2]] = area_mat[,,good_dye[2]]/area_tot[,good_dye[2]]

## Threshold abundance matrices







