## This script is used to clean up TRFLP profile data from control runs which has already been aligned to a size standard.
## It uses a peak size data table exported from PeakScanner2 to:
## 1. Remove peaks from un-used dyes and from the size standard
## 2. Remove undersized peaks due to noise (thresholding)
## 3. Remove peaks from procedural controls


options(stringsAsFactors=F)

## Define variables that will change for each analysis

# File name of PeakScanner2 txt output
run_name = 'myco_control' #'photo_control'
profiles_dir = 'myco_control_TRFLP/' #'photo_control_TRFLP/'
profiles_file = 'myco_control_TRFLP/myco_control.txt' #'photo_control_TRFLP/photo_control.txt'

# Focal dyes
good_dye = c('R','G') #c('B','Y')
names(good_dye) = c('Hpy188III','BstUI') #c('Hpy188III', 'BssKI')

# Size standard dye
size_standard_dye = 'O'
size_standard = c(as.numeric(sapply(c(0,100,200,300,400,500), function(x) x + c(14,20,40,60,80,100))),250)
size_standard = size_standard[-1]
size_standard = size_standard[order(size_standard)]

# Samples to analyze
target_samples = paste('C',1:2,'-R125-L25', sep='') # paste('P',5:8,'-R125-L50', sep='')

# Sample names for procedural controls
blank_sample = paste('L25',c('a','b','c','d'), sep='-')#'L50'
control_sample = 'PCR-R125-L25' #'PCR1-R100-L50'

# Max primer size used to eliminate peaks from primer dimers
primer_length = 20

## Read in Data

# Set data, working, and code directories
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/TRFLP/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/GitHub/myco-photo-diversity/'
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Figures/'
derived_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Derived_Data/'

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

# Remove peaks whose length is less than the length of 1.5*longest primer
peaks = subset(peaks, Size >= primer_length*1.5)

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
const_pct_mat[,size_cols] # Photobiont: Omits most size standards, Mycobiont: Omits 2 size standards
var_pct_mat[,size_cols] # Photobiont: Finds all size standards, Mycobiont: Omits 1 size standard at 500bp
use_mat[,size_cols]

# Do methods omit non-size standards?
which(const_pct_mat[,!size_cols]>0, arr.ind=T) # Photobiont & Mycobiont: Omits all non-size standard
which(var_pct_mat[,!size_cols]>0, arr.ind=T) # Photobiont: Omits all non-size standard peaks, Mycobiont: Omits all but one non-size standard peaks

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
## photobiont control: variable percentage method
## mycobiont control: variable percentage method
## substitute this method into code below


# Myco control has 4 blank samples
# Examine them 
blanks = area_mat[blank_sample, , good_dye[2]]
blanks[,colSums(blanks)>0]


## Threshold actual data
## IT IS POSSIBLE THAT THIS IS NOT REALLY REMOVING NOISE- NEED TO CHECK FINAL PROFILES AGAINST EXPECTED LENGTHS
thresh_mat1 = thresh_pct_var(area_mat[,,good_dye[1]])
thresh_mat2 = thresh_pct_var(area_mat[,,good_dye[2]])

# Try not thresholding
thresh_mat1 = area_mat[,,good_dye[1]]
thresh_mat2 = area_mat[,,good_dye[2]]

# Plot un-thresholded profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'.pdf',sep=''),
	height=1*dim(area_mat)[1], width=4)
plot_profiles(area_mat[,,good_dye[1]], blank=blank_sample, control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'.pdf',sep=''),
	height=1*dim(area_mat)[1], width=4)
plot_profiles(area_mat[,,good_dye[2]], blank=blank_sample, control=control_sample)
dev.off()

# Plot thresholded profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_threshold.pdf',sep=''),
	height=1*dim(area_mat)[1], width=4)
plot_profiles(thresh_mat1, blank=blank_sample, control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_threshold.pdf',sep=''),
	height=1*dim(area_mat)[1], width=4)
plot_profiles(thresh_mat2, blank=blank_sample, control=control_sample)
dev.off()

# After inspecting blank sample for if there is contamination that needs to be accounted for
# drop peaks less than the tallest blank peak
thresh_mat1[thresh_mat1<=max(thresh_mat1[blank_sample,])] = 0
thresh_mat2[thresh_mat2<=max(thresh_mat2[blank_sample,])] = 0

# remove blank sample from the analysis
thresh_mat1 = thresh_mat1[!(rownames(thresh_mat1) %in% blank_sample),]
thresh_mat2 = thresh_mat2[!(rownames(thresh_mat2) %in% blank_sample),]

# Drop fragment lengths that no longer contain peaks
thresh_mat1 = thresh_mat1[,colSums(thresh_mat1)>0]
thresh_mat2 = thresh_mat2[,colSums(thresh_mat2)>0]

## Align profiles using hierarchical clustering
align_mat1 = align_peaks(thresh_mat1)
align_mat2 = align_peaks(thresh_mat2)

# Plot aligned profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_aligned.pdf',sep=''),
	height=1*nrow(align_mat1), width=4)
plot_profiles(align_mat1, control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_aligned.pdf',sep=''),
	height=1*nrow(align_mat2), width=4)
plot_profiles(align_mat2, control=control_sample)
dev.off()

## Remove peaks that correspond to peaks in the procedural control
profile_mat1 = drop_peaks(align_mat1, control_sample)
profile_mat2 = drop_peaks(align_mat2, control_sample)

# Plot profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_final.pdf',sep=''),
	height=1*nrow(profile_mat1), width=4)
plot_profiles(profile_mat1)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_final.pdf',sep=''),
	height=1*nrow(profile_mat2), width=4)
plot_profiles(profile_mat2)
dev.off()

# Save profiles
write.csv(profile_mat1, paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[1],'.csv', sep=''), row.names=T)
write.csv(profile_mat2, paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[2],'.csv', sep=''), row.names=T)

## Read in cleaned profiles
profile_mat1 = read.csv(paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[1],'.csv', sep=''), row.names=1, check.names=F)
profile_mat2 = read.csv(paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[2],'.csv', sep=''), row.names=1, check.names=F)

### Compare to expected lengths
photo_lengths = read.csv(paste(derived_dir, 'expected_photobiont_strain_lengths.csv', sep=''), row.names=1)
photo_lengths = photo_lengths + 19 # Add length of ITS1 forward primer

myco_lengths = read.csv(paste(derived_dir,'expected_mycobiont_strain_lengths.csv', sep=''), row.names=1)
myco_lengths = myco_lengths + 18 # Add length of 58m2 forward primer
colnames(myco_lengths) = names(good_dye)

lengths = myco_lengths

# Define community matrix
comm = read.csv(paste(data_dir, profiles_dir,'myco_control_comm.csv', sep=''), row.names=1)

# Calculate expected lengths
len1 = lengths[rownames(comm),names(good_dye)[1]]
exp_len1 = comm * len1
len2 = lengths[rownames(comm),names(good_dye)[2]]
exp_len2 = comm * len2

# Make expected profile matrices
exp_mat1 = sapply(len1[order(len1)], function(i) 1:ncol(comm) %in% which(exp_len1 == i, arr.ind=T)[,'col'])
colnames(exp_mat1) = len1[order(len1)]
rownames(exp_mat1) = colnames(comm)
exp_mat2 = sapply(len1[order(len2)], function(i) 1:ncol(comm) %in% which(exp_len2 == i, arr.ind=T)[,'col'])
colnames(exp_mat2) = len2[order(len2)]
rownames(exp_mat2) = colnames(comm)

# Compare profiles to expected
rowSums(profile_mat1>0)
rowSums(exp_mat1)

# For photobiont and mycobiont controls, profiles are very different from expected.





