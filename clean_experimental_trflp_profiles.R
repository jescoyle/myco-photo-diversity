## This script is used to clean up TRFLP profile data from experimental runs which has already been aligned to a size standard.
## It uses a peak size data table exported from PeakScanner2 to:
## 1. Remove peaks from un-used dyes and from the size standard
## 2. Remove undersized peaks due to noise (thresholding)
## 3. Remove peaks from procedural controls
## 4. Generate a profile library from known samples


options(stringsAsFactors=F)

## Define variables that will change for each analysis

# File name of PeakScanner2 txt output
run_name = 'myco' # 'photo'
profiles_dir = 'myco_TRFLP/' # 'photo_TRFLP/'
profiles_file = paste(profiles_dir, run_name, '_trflp.txt', sep='') 

# Focal dyes
good_dye = c('R','G') #c('R','G') c('B','Y')
names(good_dye) = c('Hpy188III','BstUI') #c('Hpy188III','BstUI') c('Hpy188III', 'BssKI')

# Size standard dye
size_standard_dye = 'O'
size_standard = c(as.numeric(sapply(c(0,100,200,300,400,500), function(x) x + c(14,20,40,60,80,100))),250)
size_standard = size_standard[-1]
size_standard = size_standard[order(size_standard)]

# Sample names for procedural controls
blank_sample = 'control' # 'control' 'empty'
control_sample = 'PCR-control' # 'PCR-control' 'PCR1-control'

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

# Samples to analyze
sample_names = unique(peaks$Sample)
exp_samples = sample_names[grep('^S', sample_names)]
sp_samples = sample_names[!(sample_names %in% c(control_sample, blank_sample, exp_samples))]
target_samples = c(sp_samples, exp_samples)

# Remove peaks not assigned a size
peaks = subset(peaks, !is.na(Size))

# Remove peaks whose length is less than the minimum of the size standard
peaks = subset(peaks, Size >= min(size_standard))

# Remove peaks whose length is less than the length of 1.5*longest primer
peaks = subset(peaks, Size >= primer_length*1.5)

# Subset data to focal samples and controls
peaks = subset(peaks, Sample %in% c(sp_samples, exp_samples, blank_sample, control_sample))

## Convert dataframe to array of abundance matrices based on peak area and peak height

# Tabulate peaks
height_mat = xtabs(Height~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))
area_mat = xtabs(Area~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))

# Examine correlation between total signal and number of peaks in each sample
# This determines whether thresholding needs to account for variations in amount of DNA
par(mfrow=c(3, length(good_dye)))
par(mar=c(2,4,4,1))
for(i in good_dye){
	use_mat = area_mat[exp_samples,,i]
	npeaks = apply(use_mat, 1, function(x) sum(x>0))
	totsignal = apply(use_mat, 1, sum)
	plot(totsignal, npeaks, main=paste('Samples',i))
	mtext(paste('r =',format(cor(totsignal, npeaks),3,0, digits=3)))
}
for(i in good_dye){
	use_mat = area_mat[sp_samples,,i]
	npeaks = apply(use_mat, 1, function(x) sum(x>0))
	totsignal = apply(use_mat, 1, sum)
	plot(totsignal, npeaks, main=paste('Species',i))
	mtext(paste('r =',format(cor(totsignal, npeaks),3,0, digits=3)))
}
for(i in good_dye){
	use_mat = area_mat[c(target_samples, control_sample),,i]
	npeaks = apply(use_mat, 1, function(x) sum(x>0))
	totsignal = apply(use_mat, 1, sum)
	plot(totsignal, npeaks, main=paste('All',i))
	mtext(paste('r =',format(cor(totsignal, npeaks),3,0, digits=3)))
}

# Compare thresholding methods on size standard to be sure they don't omit a known peak but do eliminate false peaks
standard_mat = area_mat[c(target_samples,control_sample),,size_standard_dye]
std_mat_sp = area_mat[c(sp_samples,control_sample),,size_standard_dye]
std_mat_samp = area_mat[c(exp_samples,control_sample),,size_standard_dye]

use_mat = standard_mat
const_pct_mat = thresh_pct_const(use_mat)
var_pct_mat = thresh_pct_var(use_mat)

# Which columns represent size standards?
size_cols = colnames(use_mat) %in% as.character(size_standard)
colnames(use_mat)[size_cols] # There should be 35 since we dropped fragments <= 20

# Can methods find the size standards
colSums(const_pct_mat[,size_cols]>0) # Myco: Finds all, Photo: Omits 500 and 514 in most samples
colSums(var_pct_mat[,size_cols]>0) # Myco: Finds all, Photo: Finds all size standards
use_mat[,size_cols]

# Do methods omit non-size standards?
which(const_pct_mat[,!size_cols]>0, arr.ind=T) # Myco: Leaves in 5 non-size standard peaks and 1 that is in PCR control, Photo: Omits all non-size standard
which(var_pct_mat[,!size_cols]>0, arr.ind=T) # Myco: Leaves in several non-size standard peaks, Photo: Leaves in some non-size standard peaks
bad_peaks = const_pct_mat[,!size_cols][,unique(which(var_pct_mat[,!size_cols]>0, arr.ind=T)[,2])]
align_peaks(bad_peaks) 

# Compare constant and variable percentage thresholding methods
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
npeaks_stat==length(size_standard)-1 # Removes all peaks

# Examine blanks
which(area_mat[blank_sample,,good_dye[1]]>0) # Myco: 8 strong peaks
which(area_mat[blank_sample,,good_dye[2]]>0) # Myco: none

## Threshold actual data
## Do with variable percent (VP), with constant threshold (CP), with statistical method (ST)
Tmethod = 'ST'

if(Tmethod == 'CP'){
	thresh_mat1 = thresh_pct_const(area_mat[,,good_dye[1]])
	thresh_mat2 = thresh_pct_const(area_mat[,,good_dye[2]])
}
if(Tmethod == 'VP'){
	thresh_mat1 = thresh_pct_var(area_mat[,,good_dye[1]])
	thresh_mat2 = thresh_pct_var(area_mat[,,good_dye[2]])	
}
if(Tmethod == 'ST'){
	thresh_mat1 = thresh_stdev(area_mat[,,good_dye[1]], sfact=3)
	thresh_mat2 = thresh_stdev(area_mat[,,good_dye[2]], sfact=3)	
} # THIS METHOD DOES NOT WORK WELL WITH THIS DATA

# Plot un-thresholded profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+2), width=4)
plot_profiles(area_mat[c(exp_samples, blank_sample, control_sample),,good_dye[1]], blank=blank_sample, control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+2), width=4)
plot_profiles(area_mat[c(exp_samples, blank_sample, control_sample),,good_dye[2]], blank=blank_sample, control=control_sample)
dev.off()

# Plot thresholded profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_threshold_',Tmethod,'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+2), width=4)
plot_profiles(thresh_mat1[c(exp_samples, blank_sample, control_sample),], blank=blank_sample, control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_threshold_',Tmethod,'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+2), width=4)
plot_profiles(thresh_mat2[c(exp_samples, blank_sample, control_sample),], blank=blank_sample, control=control_sample)
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
align_mat1 = align_peaks(thresh_mat1, drawplot=F)
align_mat2 = align_peaks(thresh_mat2, drawplot=F)

# Plot aligned profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_aligned_',Tmethod,'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+1), width=4)
plot_profiles(align_mat1[c(exp_samples, control_sample),], control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_aligned_',Tmethod,'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+1), width=4)
plot_profiles(align_mat2[c(exp_samples, control_sample),], control=control_sample)
dev.off()

## Remove peaks that correspond to peaks in the procedural control
profile_mat1 = drop_peaks(align_mat1, control_sample) # For mycobiont: 12/206,  For photobiont: 2/77
profile_mat2 = drop_peaks(align_mat2, control_sample) # For mycobiont: 6/140, For photobiont: 5/160

# Plot profiles
pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[1],'_',Tmethod,'_final.pdf',sep=''),
	height=1*nrow(profile_mat1), width=4)
plot_profiles(profile_mat1)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_',Tmethod,'_final.pdf',sep=''),
	height=1*nrow(profile_mat2), width=4)
plot_profiles(profile_mat2)
dev.off()

# Save profiles
write.csv(profile_mat1, paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[1],'_',Tmethod,'.csv', sep=''), row.names=T)
write.csv(profile_mat2, paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[2],'_',Tmethod,'.csv', sep=''), row.names=T)


# Check whether there are any profiles without peaks
sum(rowSums(profile_mat1)==0)
sum(rowSums(profile_mat2)==0)

