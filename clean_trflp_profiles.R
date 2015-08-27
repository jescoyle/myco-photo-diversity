## This script is used to clean up TRFLP profile data which has already been aligned to a size standard.
## It uses a peak size data table exported from PeakScanner2 to:
## 1. Remove peaks from un-used dyes and from the size standard
## 2. Remove undersized peaks due to noise (thresholding)
## 3. Remove peaks from procedural controls


options(stringsAsFactors=F)

## Define variables that will change for each analysis

# File name of PeakScanner2 txt output
profiles_file = 'photo_control_TRFLP/photo_control.txt'

# Focal dyes
good_dye = c('B','Y')
names(good_dye) = c('Hpy188III', 'BSSKI')

# Size standard dye
size_standard_dye = 'O'
size_standard = c(as.numeric(sapply(c(0,100,200,300,400,500), function(x) x + c(14,20,40,60,80,100))),250)
size_standard = size_standard[-1]
size_standard = size_standard[order(size_standard)]

# Samples to analyze
target_samples = paste('P',5:8,'-R125-L50', sep='')

# Sample names for procedural controls
blank_sample = 'L50'
control_sample = 'PCR1-R100-L50'

## Read in Data

# Set data, working, and code directories
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/TRFLP/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/GitHub/myco-photo-diversity/'

setwd(working_dir)

# Read in functions for the analysis of TRFLP data
source(paste(git_dir,'TRFLP_functions.R', sep=''))

# Read in PeakScanner2 file and re-name columns
peaksc_raw = read.table(paste(data_dir,profiles_file, sep=''), sep='\t', header=T, check.names=F)
peaks = peaksc_raw[,c('Sample Name','Dye/Sample Peak','Size','Height','Width in BP', 'Area in BP','UD3')]
names(peaks) = c('Sample','DyeNum','Size','Height','Width', 'Area','strain')

# Split Dye and Peak numbers
peaks$Dye = substr(peaks$DyeNum, 1, 1)
peaks$Num = as.numeric(extract_vals(peaks$DyeNum, regexpr('[A-Z], ([0-9]+)', peaks$DyeNum, perl=T)))

# Remove peaks not assigned a size
peaks = subset(peaks, !is.na(Size))

# Remove peaks whose length is less than the minimum of the size standard
peaks = subset(peaks, Size >= min(size_standard))

# Subset data to focal samples and controls
peaks = subset(peaks, Sample %in% c(target_samples, blank_sample, control_sample))

## Convert dataframe to array of abundance matrices based on peak area and peak height

# Tabulate peaks
height_mat = xtabs(Height~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))
area_mat = xtabs(Area~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))

# Examine correlation between total signal and number of peaks in each sample
# This determines whether thresholding needs to account for variations in amount of DNA
par(mfrow=c(1, length(good_dye)))
for(i in good_dye){
	use_mat = area_mat[target_samples,,i]
	npeaks = apply(use_mat, 1, function(x) sum(x>0))
	totsignal = apply(use_mat, 1, sum)
	plot(totsignal, npeaks, main=i)
	mtext(paste('r =',format(cor(totsignal, npeaks),3,0, digits=3)))
}

# Compare thresholding methods on size standard to be sure they don't omit a known peak but do eliminate false peaks
standard_mat = area_mat[c(target_samples,control_sample),,size_standard_dye]

use_mat = standard_mat
const_pct_mat = thresh_pct_const(use_mat)
var_pct_mat = thresh_pct_var(use_mat)

# Which columns represent size standards?
size_cols = colnames(use_mat) %in% as.character(size_standard)

# Can methods find the size standards
const_pct_mat[,size_cols] # Omits most size standards
var_pct_mat[,size_cols] # Finds all size standards

# Do methods omit non-size standards?
which(const_pct_mat[,!size_cols]>0, arr.ind=T) # Omits all non-size standard
which(var_pct_mat[,!size_cols]>0, arr.ind=T) # Omits all non-size standard peaks

# Compare constant and variable percentage thresholding methods
const_pct_mat = thresh_pct_const(use_mat)
var_pct_mat = thresh_pct_var(use_mat)
par(mfrow=c(1,3))
npeaks = apply(use_mat, 1, function(x) sum(x>0))
totsignal = apply(use_mat, 1, sum)
plot(totsignal, npeaks, main='Unthresholded Data')
npeaks = apply(const_pct_mat, 1, function(x) sum(x>0))
totsignal = apply(const_pct_mat, 1, sum)
plot(totsignal, npeaks, main='Constant Percentage')
npeaks = apply(var_pct_mat, 1, function(x) sum(x>0))
totsignal = apply(var_pct_mat, 1, sum)
plot(totsignal, npeaks, main='Variable Percentage')

# Statistical thresholding method does not work well with this data.
use_mat = standard_mat
npeaks_stat = sapply(seq(3, 21,2), function(x){
	thresh_mat = thresh_stdev(use_mat,x)
	apply(thresh_mat, 1, function(x) sum(x>0))
})
npeaks_stat==length(size_standard) # Never removes all non-size standard peaks from PCR run

## Decide which threshold method to use:
## photobiont control: variable threshold method
## substitute this method into code below


## Threshold actual data
thresh_mat1 = thresh_pct_var(area_mat[,,good_dye[1]])
thresh_mat2 = thresh_pct_var(area_mat[,,good_dye[2]])


## Standardize abundance matrices by total area or total height
height_tot = apply(height_mat, c(1,3), sum)
height_mat[,,good_dye[1]] = height_mat[,,good_dye[1]]/height_tot[,good_dye[1]]
height_mat[,,good_dye[2]] = height_mat[,,good_dye[2]]/height_tot[,good_dye[2]]

area_tot = apply(area_mat, c(1,3), sum)
area_mat[,,good_dye[1]] = area_mat[,,good_dye[1]]/area_tot[,good_dye[1]]
area_mat[,,good_dye[2]] = area_mat[,,good_dye[2]]/area_tot[,good_dye[2]]

## Threshold abundance matrices









