## This script is used to clean up TRFLP profile data from experimental runs which has already been aligned to a size standard.
## It uses a peak size data table exported from PeakScanner2 to:
## 1. Remove peaks from un-used dyes and from the size standard
## 2. Remove undersized peaks due to noise (thresholding)
## 3. Remove peaks from procedural controls
## 4. Generate a profile library from known samples

options(stringsAsFactors=F)

## Define variables that will change for each analysis

# File name of PeakScanner2 txt output
run_name = 'photo_Dec2015' #'myco' # 'photo'  'photo_Dec2015'
profiles_dir = 'photo_Dec2015/' #'photo_Dec2015/' # 'photo_TRFLP/'
profiles_file = paste0(profiles_dir, run_name, c('_run1.txt', '_run2.txt')) # paste0(profiles_dir, run_name, c('_run1.txt', '_run2.txt')) # paste0(profiles_dir, run_name, '_trflp.txt') 

# Focal dyes
good_dye = 'Y' #'Y' c('R','G') c('B','Y')
names(good_dye) = 'MspI' #'MspI' c('Hpy188III','BstUI') c('Hpy188III', 'BssKI')

# Size standard dye
size_standard_dye = 'O'
size_standard = c(as.numeric(sapply(c(0,100,200,300,400,500), function(x) x + c(14,20,40,60,80,100))),250)
size_standard = size_standard[-1]
size_standard = size_standard[order(size_standard)]

# Sample names for procedural controls
blank_sample = 'empty' # 'control' 'empty'
control_sample = c('PCR1-control', 'PCR2-control') # c('PCR1-control', 'PCR2-control') 'PCR-control' 'PCR1-control'

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

###########################################################
### Assess repeatability using duplicate runs from Dec 2015

align_list = sapply(profiles_file, function(f){
	peaks = read_peakscanner(paste0(data_dir, f), size_standard, primer_length)
	area_mat = xtabs(Area~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))
	align_mat = align_peaks(area_mat[,,good_dye], F)
	align_mat
}, simplify='array')

# Combine into array
frag_lens = unique(c(dimnames(align_list[[1]])[[2]], dimnames(align_list[[2]])[[2]]))
frag_lens = frag_lens[order(as.numeric(frag_lens))]
align_mat = sapply(frag_lens, function(x){
	col1 = if(x %in% dimnames(align_list[[1]])[[2]]){ align_list[[1]][,x]} else { rep(0, nrow(align_list[[1]]))}
	col2 = if(x %in% dimnames(align_list[[2]])[[2]]){ align_list[[2]][,x]} else { rep(0, nrow(align_list[[1]]))}
	cbind(col1, col2)
}, simplify='array')


# Image plot of whether two runs give the same or different peaks

presence = apply(align_mat, c(1,3), function(x){
	x = as.numeric(x > 0)
	if(sum(x)==0) x = NA
	sum(x*c(-1,1))
})

image(t(presence), col=colorRampPalette(c('blue','green','yellow'))(3))#, col = use_col[factor(presence)])

std_mat = apply(align_mat, c(1,2), function(x) x/sum(x)) 
std_mat[is.na(std_mat)] = 0

color_mat = apply(std_mat, c(1,2), function(x) rgb(1-x[1],1,1-x[2]))
#color_mat = apply(std_mat, c(1,2), function(x) rgb(x[1],0,x[2]))
nsamps = ncol(color_mat)
npeaks = nrow(color_mat)

png(paste0(fig_dir, 'Dec2015_peak_replication.png'), height=2000, width=3000)
plot(1,1, type='n', axes=F, xlab='Peaks', ylab='Samples',xlim=c(-1, npeaks+1), ylim=c(-1, nsamps+1))
for(j in 1:nsamps){
for(i in 1:npeaks){
	rect(i-1, j-1, i, j, col=color_mat[i,j], border=NA)
}}
axis(2, at=(1:nsamps)-0.5, labels=colnames(color_mat), las=1, tick=F)
dev.off()

# This image plot makes it pretty clear that peaks within 1 unit should probably be combined
# Try re-aligning with bin size of 2

# Align peaks from replicate runs to each other
peaks_list = sapply(profiles_file, function(f){
	peaks = read_peakscanner(paste0(data_dir, f), size_standard, primer_length)
	area_mat = xtabs(Area~Sample+Size+Dye, data=subset(peaks, Dye %in% c(good_dye,size_standard_dye)))
	area_mat
}, simplify='array')

combined_mat = cbind(peaks_list[[1]][,,good_dye], peaks_list[[2]][,,good_dye])
sizes = as.numeric(colnames(combined_mat))
dmat = dist(sizes)
cl = hclust(dmat, method='average')
bins = cutree(cl, h=2) # Determine bin width here
bins1 = bins[1:ncol(peaks_list[[1]][,,good_dye])]
bins2 = bins[(ncol(peaks_list[[1]][,,good_dye])+1):length(bins)]

align_mat = array(0, dim=c(dim(peaks_list[[1]])[1], max(bins), 2))
for(i in 1:max(bins)){
	if(i %in% bins1){
		merge_cols = which(bins1==i)
		if(length(merge_cols)>1){
			align_mat[,i,1] = rowSums(peaks_list[[1]][, merge_cols, good_dye])
		} else {
			align_mat[,i,1] = peaks_list[[1]][, merge_cols, good_dye]
		}
	}
	if(i %in% bins2){
		merge_cols = which(bins2==i)
		if(length(merge_cols)>1){
			align_mat[,i,2] = rowSums(peaks_list[[2]][, merge_cols, good_dye])
		} else {
			align_mat[,i,2] = peaks_list[[2]][, merge_cols, good_dye]
		}
	}
}

# Use mean of sizes to label bins
bin_names = sapply(1:max(bins), function(x) round(mean(sizes[bins==x]),0))
dimnames(align_mat) = list(dimnames(peaks_list[[1]])$Sample, bin_names, 1:2)

# For each run, drop peaks that are smaller than the largest peak in the blank sample
clean_mat = align_mat
clean_mat[,,1][clean_mat[,,1]<=max(align_mat[blank_sample,,1])] = 0
clean_mat[,,2][clean_mat[,,2]<=max(align_mat[blank_sample,,2])] = 0

# Remove peaks that show up in both runs from controls
contaminated = apply(clean_mat[control_sample,,]>0, 1:2, sum)
drop_peaks = apply(contaminated==2, 2, function(x) sum(x) > 0)
clean_mat[,drop_peaks,] = 0

# Standardize by total area
std_mat = apply(clean_mat, c(1,3), function(x) x/sum(x)) 
std_mat[is.na(std_mat)] = 0

color_mat = apply(std_mat, c(1,2), function(x) rgb(1-x[1],1,1-x[2]))
nsamps = ncol(color_mat)
npeaks = nrow(color_mat)

png(paste0(fig_dir, 'Dec2015_peak_replication_bin2.png'), height=3500, width=4500)
plot(1,1, type='n', axes=F, xlab='Peaks', ylab='Samples',xlim=c(-1, npeaks+1), ylim=c(-1, nsamps+1))
for(j in 1:nsamps){
for(i in 1:npeaks){
	rect(i-1, j-1, i, j, col=color_mat[i,j], border=NA)
}}
abline(h=0:nsamps, col='gray')
abline(v=0:npeaks, col='gray')
axis(2, at=(1:nsamps)-0.5, labels=colnames(color_mat), las=1, tick=F, cex=2, pos=0)
axis(1, at=(1:npeaks)-0.5, labels=rownames(color_mat), las=2, tick=F, cex=2, pos=0)
axis(3, at=(1:npeaks)-0.5, labels=rownames(color_mat), las=2, tick=F, cex=2, pos=nsamps)
axis(4, at=(1:nsamps)-0.5, labels=colnames(color_mat), las=1, tick=F, cex=2, pos=npeaks)
dev.off()

# Presence-only
both_present = apply(align_mat, 1:2, function(x) sum(x>0)==2)
png(paste0(fig_dir, 'Dec2015_peak_replication_bin2_presence.png'), height=3500, width=4500)
plot(1,1, type='n', axes=F, xlab='Peaks', ylab='Samples',xlim=c(-1, npeaks+1), ylim=c(-1, nsamps+1))
for(j in 1:nsamps){
for(i in 1:npeaks){
	rect(i-1, j-1, i, j, col=as.numeric(both_present[j,i]), border=NA)
}}
abline(h=0:nsamps, col='gray')
abline(v=0:npeaks, col='gray')
axis(2, at=(1:nsamps)-0.5, labels=colnames(color_mat), las=1, tick=F, cex=2, pos=0)
axis(1, at=(1:npeaks)-0.5, labels=rownames(color_mat), las=2, tick=F, cex=2, pos=0)
axis(3, at=(1:npeaks)-0.5, labels=rownames(color_mat), las=2, tick=F, cex=2, pos=nsamps)
axis(4, at=(1:nsamps)-0.5, labels=colnames(color_mat), las=1, tick=F, cex=2, pos=npeaks)
dev.off()

# For each sample, what is the largest peak not replicated across runs?
apply(std_mat, 2, function(x){
	one_absent = rowSums(x > 0)==1
	max(x[one_absent])
})
# For the most part all major peaks are observed except: cat_nig, fla_cap, S29, S47, S52

# Plot profiles from both runs
exp_samples = dimnames(std_mat)[[2]][!(dimnames(std_mat)[[2]] %in% c(blank_sample, control_sample))]

plot_data = std_mat[,exp_samples,]
N = dim(plot_data)[2]

pdf(paste0(fig_dir, 'Dec2015_peak_replication_profiles_bin2.pdf'), height=12, width=5)
par(mfrow=c(5,1))
par(mar=c(2,4,0.5,1))
par(lend=1)
for(i in 1:N){
	this_pf = plot_data[,i,]
	this_pf[this_pf==0] = NA
	plot(as.numeric(rownames(this_pf)), this_pf[,1], type='h', , col=rgb(0,0,1),
		ylim=c(-1,1), ylab=dimnames(plot_data)[[2]][i], xlab='', las=1, axes=F, bg=1)
	points(as.numeric(rownames(this_pf)), -this_pf[,2], type='h', col=rgb(1,0,0))
	axis(2, las=1)
	if((i%%5 == 0)|(i==N)) axis(1)
}
dev.off()

# Re-plot just the species with expected fragment lengths
sp_samples = exp_samples[grep('[a-z]+-[a-z]+', exp_samples)]

# Read in expected lengths
photo_lengths = read.csv(paste0(derived_dir, 'expected_photobiont_strain_length_full.csv'), row.names=1)
strain_names = read.csv(paste0(derived_dir, 'photobiont_strain_names.csv'))

pdf(paste0(fig_dir, 'Dec2015_peak_replication_species_profiles_bin2.pdf'), height=12, width=5)
par(mfrow=c(5,1))
par(mar=c(2,4,0.5,1))
par(lend=1)
for(i in 1:length(sp_samples)){
	this_sp = sp_samples[i]
	this_pf = plot_data[,this_sp,]
	this_pf[this_pf==0] = NA
	plot(as.numeric(rownames(this_pf)), this_pf[,1], type='h', , col=rgb(0,0,1),
		ylim=c(-1,1), ylab=this_sp, xlab='', las=1, las=1, bg=1)
	points(as.numeric(rownames(this_pf)), -this_pf[,2], type='h', col=rgb(1,0,0))
	this_strain = unique(subset(strain_names, SpID == this_sp)$Strain_photo)
	if(length(this_strain)>0) points(photo_lengths[this_strain,'MspI'], rep(0, length(this_strain)), pch=1)
}
dev.off()


## Threshold by excluding peaks not present in both runs
both_present = apply(clean_mat, 1:2, function(x) sum(x>0)==2)
thresh_mat = apply(clean_mat, 1:2, function(x){
	if(x[1]>0 & x[2]>0){
		x
	} else {
		c(0,0)
	}	
})

# Save profiles from runs separately
write.csv(thresh_mat[1,,], paste0(derived_dir, 'profiles_photo_MspI_compareruns-run1.csv', sep=''), row.names=T)
write.csv(thresh_mat[2,,], paste0(derived_dir, 'profiles_photo_MspI_compareruns-run2.csv', sep=''), row.names=T)

# Save combined profiles
combined = apply(thresh_mat, 2:3, sum)
write.csv(combined, paste0(derived_dir, 'profiles_photo_MspI_compareruns-combined.csv', sep=''), row.names=T)


## Threshold each run separately using methods below that were use in previous runs
## Used code below:

######################################################################
### Use this code when analyzing a single run
# Read in PeakScanner2 file and re-name columns

#run_name = 'photo_Dec2015-1'
peaks = read_peakscanner(paste0(data_dir, profiles_file[2]), size_standard, primer_length)

# Samples to analyze
sample_names = unique(peaks$Sample)
exp_samples = sample_names[grep('^S', sample_names)]
sp_samples = sample_names[!(sample_names %in% c(control_sample, blank_sample, exp_samples))]
target_samples = c(sp_samples, exp_samples)

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
colSums(const_pct_mat[,size_cols]>0) # Myco: Finds all, Photo: Omits 500 and 514 in most samples, Photo-Dec2015: Omits several in all samples
colSums(var_pct_mat[,size_cols]>0) # Myco: Finds all, Photo: Finds all size standards, Photo-Dec2015: Finds all
use_mat[,size_cols]

# Do methods omit non-size standards?
which(const_pct_mat[,!size_cols]>0, arr.ind=T) # Myco: Leaves in 5 non-size standard peaks and 1 that is in PCR control, Photo: Omits all non-size standard, Photo-Dec2015: Omits all non-size standard
which(var_pct_mat[,!size_cols]>0, arr.ind=T) # Myco: Leaves in several non-size standard peaks, Photo: Leaves in some non-size standard peaks, Photo-Dec2015: Leaves in several non-size standard peaks
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
area_mat[blank_sample,which(area_mat[blank_sample,,good_dye[1]]>0),good_dye[1]] # Myco: 8 strong peaks, Photo: none, Photo Dec2015-1: 7 small, Dec2015-2: 3 small
area_mat[blank_sample,which(area_mat[blank_sample,,good_dye[2]]>0),good_dye[2]] # Myco: none, Photo: one small peak

# After inspecting blank sample for if there is contamination that needs to be accounted for
# drop peaks less than the tallest blank peak
clean_mat1 = area_mat[,,good_dye[1]]
clean_mat2 = area_mat[,,good_dye[2]]
clean_mat1[clean_mat1<=max(clean_mat1[blank_sample,])] = 0
clean_mat2[clean_mat2<=max(clean_mat2[blank_sample,])] = 0

# remove blank sample from the analysis
clean_mat1 = clean_mat1[!(rownames(clean_mat1) %in% blank_sample),]
clean_mat2 = clean_mat2[!(rownames(clean_mat2) %in% blank_sample),]


## Threshold actual data
## Do with variable percent (VP), with constant threshold (CP), with statistical method (ST)
Tmethod = 'CP'

if(Tmethod == 'CP'){
	thresh_mat1 = thresh_pct_const(clean_mat1)
	thresh_mat2 = thresh_pct_const(clean_mat2)
}
if(Tmethod == 'VP'){
	thresh_mat1 = thresh_pct_var(clean_mat1)
	thresh_mat2 = thresh_pct_var(clean_mat2)	
}
if(Tmethod == 'ST'){
	thresh_mat1 = thresh_stdev(clean_mat1, sfact=3)
	thresh_mat2 = thresh_stdev(clean_mat2, sfact=3)	
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
plot_profiles(thresh_mat1[c(exp_samples, control_sample),], control=control_sample)
dev.off()

pdf(paste(fig_dir,'profiles_',run_name,'_',names(good_dye)[2],'_threshold_',Tmethod,'_samples.pdf',sep=''),
	height=1*(length(exp_samples)+2), width=4)
plot_profiles(thresh_mat2[c(exp_samples, control_sample),], control=control_sample)
dev.off()

# Drop fragment lengths that no longer contain peaks
thresh_mat1 = thresh_mat1[,colSums(thresh_mat1)>0]
thresh_mat2 = thresh_mat2[,colSums(thresh_mat2)>0]

## Align profiles using hierarchical clustering
# Using bin size 2 to conform with MspI runs in Dec 2015
align_mat1 = align_peaks(thresh_mat1, drawplot=F, bin_size=2)
align_mat2 = align_peaks(thresh_mat2, drawplot=F, bin_size=2)

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
profile_mat1 = drop_peaks(align_mat1, control_sample) # For mycobiont: 12/206,  For photobiont: 2/74, For photobiont Dec2015-1: 4/110 Dec2015-2: 5/103
profile_mat2 = drop_peaks(align_mat2, control_sample) # For mycobiont: 6/140, For photobiont: 5/121

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

