## This script is used to analyze photobiont T-RFLP profiles generated from the Photobiont-Mycobiont Diversity Project

## TO DO: COMPARE RICHNESS AND COMMUNITY-ENV CORRELATIONS FROM DIFFERENT RESTRICTION ENZYMES AND THRESHOLDING METHODS


options(stringsAsFactors=F)

library(reshape)
library(cluster)

# Set data, working, and code directories
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/TRFLP/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/GitHub/myco-photo-diversity/'
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Figures/'
derived_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Derived_Data/'

setwd(working_dir)

# Load functions
source(paste(git_dir,'TRFLP_functions.R', sep=''))

# Set working colors
mycolor = read.csv('../blue2red_10colramp.txt')[9:2,]
mycolor = rgb(mycolor, maxColorValue=255)

## Define variables that will change for each analysis

# File name of PeakScanner2 txt output
run_name =  'photo' # 'photo' 'myco'
profiles_dir = 'photo_TRFLP/' # 'photo_TRFLP/' 'myco_TRFLP/'
profiles_file = paste(profiles_dir, run_name, '_trflp.txt', sep='') 

# Focal dyes
good_dye = c('B','Y') #c('R','G') c('B','Y')
names(good_dye) = c('Hpy188III', 'BssKI') #c('Hpy188III','BstUI') c('Hpy188III', 'BssKI')

# Thresholding method
Tmethod = 'CP'#'VP'

# Read in expected lengths
#photo_lengths = read.csv(paste(derived_dir, 'expected_photobiont_strain_lengths.csv', sep=''), row.names=1)
#photo_lengths = photo_lengths + 19 + 20# Add length of ITS1 + ITS2 primers
photo_lengths = read.csv(paste(derived_dir, 'expected_photobiont_strain_lengths_manual.csv', sep=''), row.names=1)

myco_lengths = read.csv(paste(derived_dir,'expected_mycobiont_strain_lengths.csv', sep=''), row.names=1)
# MYCO: DO NOT ADD PRIMER LENGTHS TO EXPECTED LENGTHS BECAUSE SEQS INCLUDE PRIMERS.
colnames(myco_lengths) = names(good_dye)

# Read in species that strains are associated with
myco_names = read.csv(paste(derived_dir, 'mycobiont_strain_names.csv', sep=''))
rownames(myco_names) = myco_names$strainID
photo_names = read.csv(paste(derived_dir, 'photobiont_strain_names.csv', sep=''))

photo_lengths$Strain_photo = rownames(photo_lengths)
photo_lengths = merge(photo_lengths, photo_names, all.x=T, all.y=F)

## Read in profile matrices
profile_mat1 = read.csv(paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[1],'_',Tmethod,'.csv', sep=''), row.names=1, check.names=F)
profile_mat2 = read.csv(paste(derived_dir, 'profiles_',run_name,'_',names(good_dye)[2],'_',Tmethod,'.csv', sep=''), row.names=1, check.names=F)

# Samples to analyze
sample_names = rownames(profile_mat1)
exp_samples = sample_names[grep('^S', sample_names)]
sp_samples = sample_names[!(sample_names %in% exp_samples)]


### Check whether experimental samples match known species profiles ###
sp_mat1 = profile_mat1[sp_samples,]
sp_mat2 = profile_mat2[sp_samples,]

# Calculate standardized peaks
sp_mat1_std = sp_mat1 / rowSums(sp_mat1)
sp_mat2_std = sp_mat2 / rowSums(sp_mat2)

# Find tallest peak in each species sample and save in library
find_largest_peak = function(x) which(x==max(x))

sp_library = data.frame(sp_samples, 
	colnames(sp_mat1)[apply(sp_mat1, 1, find_largest_peak)],
	colnames(sp_mat1)[apply(sp_mat2, 1, find_largest_peak)])
names(sp_library) = c('SpID',names(good_dye))

# Plot ranks of peaks
pdf(paste(fig_dir, run_name, '_species_peak_ranks_', names(good_dye)[1], '_',Tmethod,'_.pdf', sep=''), height=8, width=6)
par(mfrow=c(4,3))
par(mar=c(2,2,1,1))
for(s in sp_samples){
	these_peaks = sp_mat1_std[s,]
	these_peaks = these_peaks[these_peaks>0]
	plot(rank(these_peaks), these_peaks)
	ordered_peaks = these_peaks[order(these_peaks, decreasing=T)]
	cum_pct = sapply(1:length(ordered_peaks), function(i) sum(ordered_peaks[1:i]))
	cutoff = which(cum_pct==min(cum_pct[cum_pct>.5]))	
	abline(v=length(these_peaks)-cutoff+1, col=2)
	cutoff = which(cum_pct==min(cum_pct[cum_pct>.9]))	
	abline(v=length(these_peaks)-cutoff+1, col=4)
	mtext(s, 3, -2, adj=0.1)
}
dev.off()

# List all peaks and find which peaks within top 50% of cumulative area
sp_peaks1 = lapply(sp_samples, function(s){
	these_peaks = sp_mat1_std[s,]
	these_sizes = names(these_peaks)[these_peaks>0]
	these_peaks = these_peaks[these_peaks>0]
	ordered_peaks = these_peaks[order(these_peaks, decreasing=T)]
	ordered_sizes = these_sizes[order(these_peaks, decreasing=T)]
	cum_pct = sapply(1:length(ordered_peaks), function(i) sum(ordered_peaks[1:i]))
	cutoff = which(cum_pct==min(cum_pct[cum_pct>.5]))
	list(peaks=data.frame(bp=ordered_sizes, size=ordered_peaks), cut50=cutoff)
})
names(sp_peaks1)=sp_samples

sp_peaks2 = lapply(sp_samples, function(s){
	these_peaks = sp_mat2_std[s,]
	these_sizes = names(these_peaks)[these_peaks>0]
	these_peaks = these_peaks[these_peaks>0]
	ordered_peaks = these_peaks[order(these_peaks, decreasing=T)]
	ordered_sizes = these_sizes[order(these_peaks, decreasing=T)]
	cum_pct = sapply(1:length(ordered_peaks), function(i) sum(ordered_peaks[1:i]))
	cutoff = which(cum_pct==min(cum_pct[cum_pct>.5]))
	list(peaks=data.frame(bp=ordered_sizes, size=ordered_peaks), cut50=cutoff)
})
names(sp_peaks2)=sp_samples

# Compare library to expected sequence lengths
lengths = photo_lengths # myco_lengths
#rownames(lengths) = myco_names[rownames(lengths),'SpID']

sp_exp = lengths[sp_library$SpID,]

plot.new()
plot.window(xlim=c(0,600), ylim=c(0,600), xlab=names(good_dye)[1], ylab=names(good_dye)[2])
arrows(sp_exp[,1], sp_exp[,2], as.numeric(sp_library[,names(good_dye)[1]]), as.numeric(sp_library[,names(good_dye)[2]]))
axis(1, las=1)
axis(2, las=1)
mtext(names(good_dye)[1], 1, 2.5)
mtext(names(good_dye)[2], 2, 3)

# Find rank of expected length in actual species profiles
# Only works for mycobionts
# 12/6/15: May need to go one-by-one for photobionts because multiple matches for lec_str and look for partial digest
matched_peaks = array(NA, dim=c(2,length(sp_samples), 6),
	dimnames = list(Dye=names(good_dye), SpID=sp_samples, Statistic=c('bp_exp','bp_found','size','bp_diff','rank','perc')))

for(s in sp_samples){
	if(s %in% rownames(lengths)| s %in% lengths$SpID){
	peaks1 = sp_peaks1[[s]]$peaks
	peaks2 = sp_peaks2[[s]]$peaks
	if(run_name=='myco') exp_lens = lengths[s,]
	if(run_name=='photo'){
		exp_lens = subset(lengths, SpID==s)[,names(good_dye)]
	}	

	diffs1 = abs(as.numeric(peaks1$bp) - exp_lens[,1])
	rank1 = which(diffs1 == min(diffs1))
	if(length(rank1)>1) rank1 = rank1[1]
	perc1 = sum(peaks1[1:rank1,'size'])
	matched_peaks[names(good_dye)[1], s, ] = c(exp_lens[,1], as.numeric(peaks1[rank1,]), diffs1[rank1], rank1, perc1)
	
	diffs2 = abs(as.numeric(peaks2$bp) - exp_lens[,2])
	rank2 = which(diffs2 == min(diffs2))
	if(length(rank2)>1) rank2 = rank2[1]
	perc2 = sum(peaks2[1:rank2,'size'])
	matched_peaks[names(good_dye)[2], s, ] = c(exp_lens[,2], as.numeric(peaks2[rank2,]), diffs1[rank2], rank2, perc2)

	#print(sp_peaks1[[s]]$peaks)
	#print(sp_peaks2[[s]]$peaks)
	#print(lengths[s,])
	}
}

matched_peaks_df=cast(melt(matched_peaks), Dye+SpID~Statistic)

write.csv(matched_peaks_df, paste(fig_dir, run_name, '_match_peaks_',Tmethod,'.csv', sep=''), row.names=F)

# Examine species one-by-one
i = 40
sp = sp_samples[i]
subset(lengths, SpID==sp)
sp_peaks1[[sp]]
sp_peaks2[[sp]]

# Match species peaks to sample peaks to see whether species are in samples or match other species
match_array1 = array(NA, dim=c(length(sp_samples), nrow(profile_mat1), 3), 
	dimnames=list(SpID=sp_samples, SampID=rownames(profile_mat1), Statistic=c('Pct_match','N_match','N_peaks')))
for(i in sp_samples){
for(j in rownames(profile_mat1)){
	match_array1[i,j,] = match_peaks(profile_mat1[i,], profile_mat1[j,])
}}

match_array2 = array(NA, dim=c(length(sp_samples), nrow(profile_mat2), 3), 
	dimnames=list(SpID=sp_samples, SampID=rownames(profile_mat2), Statistic=c('Pct_match','N_match','N_peaks')))
for(i in sp_samples){
for(j in rownames(profile_mat2)){
	match_array2[i,j,] = match_peaks(profile_mat2[i,], profile_mat2[j,])
}}

# Save arrays
save(match_array1, match_array2, file=paste(derived_dir,run_name,'_matched_profile_arrays_',Tmethod,'.RData', sep=''))
load(paste(derived_dir,run_name,'_matched_profile_arrays_',Tmethod,'.RData', sep=''))

# Examine which peaks species have in common
image(match_array1[sp_samples,sp_samples,'Pct_match'])
sp_dist = dist(t(match_array1[sp_samples,sp_samples,'Pct_match']))
sp_clust = hclust(sp_dist, method='average')
plot(sp_clust)

## Compare known community composition to species profiles
# Note: 'lec-sp2' was IDed as 'lec-str'.
# Mycobionts: have 'lec-sp2' only so change name of community matrix to 'lec-sp2'
# Photobionts: have 'lec-str' only so leave name in community matrix alone
# Photobionts: have 'usn-sp2' which should be 'usn-str' in community matrix, but these have different photobiont strains

# Load community composition data and convert to sampXsp matrix
comm_df = read.csv('../Canopy Functional Traits/Data/SQLite Tables/community.csv')
sampXsp = 1*(xtabs(~SampID+TaxonID, data=comm_df)>0)

# Convert column names and subset to species for which we have profiles
colnames(sampXsp) = tolower(sub('_', '-', colnames(sampXsp)))
rownames(sampXsp) = paste('S', rownames(sampXsp), sep='')
if(run_name=='myco') colnames(sampXsp)[colnames(sampXsp)=='lec-str'] = 'lec-sp2'

# Add in S61 + S62 where both samples were combined
sampXsp = rbind(sampXsp, 1*(colSums(sampXsp[c('S61','S62'),])>0))
rownames(sampXsp)[nrow(sampXsp)] = 'S61-S62'

# Subset
comm_sp = sp_samples[sp_samples %in% colnames(sampXsp)]
sampXsp_sub = sampXsp[exp_samples,comm_sp] # Photobionts

# Find difference between proportion of a species' profile that is in a sample and whether the species is in the sample (0,1)
# Will be 0 when full profile in sample and species in sample
# Will be positive when species in sample and some of profile in sample
# Will be negative when species not in sample but some of profile is in sample
# Magnitude gives dissimilarity (0, 1)

spsamp_dist1 = sampXsp_sub - t(match_array1[comm_sp,exp_samples,'Pct_match']) 
spsamp_dist2 = sampXsp_sub - t(match_array2[comm_sp,exp_samples,'Pct_match']) 

Nsamp = nrow(sampXsp_sub)
Nsp = ncol(sampXsp_sub)

use_col = c(mycolor[1:4], '#FFFFFF', mycolor[5:8])


pdf(paste(fig_dir, run_name, '_compare_sp_samp_profiles_',Tmethod,'.pdf', sep=''), height=12, width=20)
layout(matrix(1:8, nrow=2, byrow=T), widths=c(0.3,.3,.3,.1))
par(mar=c(5,5,2,1))

# Presence-Absence community matrix
image(1:Nsp, 1:Nsamp, t(sampXsp_sub), axes=F, xlab='', ylab='', col=c('white','black'), main='Presence/Absence')
axis(1, at=1:Nsp, labels=comm_sp, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Composition based on profiles from Dye1
image(1:Nsp, 1:Nsamp, match_array1[comm_sp,exp_samples,'Pct_match'], axes=F, xlab='', ylab='', main=names(good_dye)[1],
	col=colorRampPalette(c('white','black'))(10), breaks=seq(0,1,.1))
axis(1, at=1:Nsp, labels=comm_sp, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Composition based on profiles from Dye2
image(1:Nsp, 1:Nsamp, match_array2[comm_sp,exp_samples,'Pct_match'], axes=F, xlab='', ylab='', main=names(good_dye)[2],
	col=colorRampPalette(c('white','black'))(10), breaks=seq(0,1,.1))
axis(1, at=1:Nsp, labels=comm_sp, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# BW legend
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
plotColorRamp(c('white','black'), 10, c(.1,.2,.3,.9), labels=seq(0,1,.1))


# Leave a blank plot
par(mar=c(5,5,2,1))
plot.new()

# Dissimilarity based on profile from Dye1
image(1:Nsp, 1:Nsamp, t(spsamp_dist1), axes=F, xlab='', ylab='', 
	col=colorRampPalette(use_col)(20), breaks=seq(-1,1,.1))
axis(1, at=1:Nsp, labels=comm_sp, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Dissimilarity based on profile from Dye2
image(1:Nsp, 1:Nsamp, t(spsamp_dist2), axes=F, xlab='', ylab='', 
	col=colorRampPalette(use_col)(20), breaks=seq(-1,1,.1))
axis(1, at=1:Nsp, labels=comm_sp, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Color legend
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
plotColorRamp(use_col, 20, c(.1,.2,.3,.9), labels=seq(-1,1,.2))
dev.off()


### Analyze community composition based on profile peaks ###

library(vegan)

# Standardize to relative peak area
mat1_std = as.matrix(profile_mat1/rowSums(profile_mat1))
mat2_std = as.matrix(profile_mat1/rowSums(profile_mat2))

# Hellinger transformation
mat1_hel = sqrt(mat1_std)
mat2_hel = sqrt(mat2_std)

# Compare distance matrices between two dyes
dmat1 = as.matrix(dist(mat1_hel))
dmat2 = as.matrix(dist(mat2_hel))
par(mfrow=c(1,2))
image(dmat1)
image(dmat2)

# Unconstrained ordination (PCA on Hellinger transformed data)
rda0_1 = rda(mat1_hel~1)
rda0_2 = rda(mat2_hel~1)

# Define matrix of explanatory variables
samples = read.csv(paste(derived_dir, 'samples.csv', sep=''))
loggers = read.csv(paste(derived_dir, 'loggerdata.csv', sep=''))
samples = merge(samples, loggers)
rownames(samples) = paste('S', samples$SampID, sep='')

# Calculate richness from community data matrix
rich = rowSums(sampXsp>0)
samples$S = rich[rownames(samples)]

# Subset env data to samples and variables of interest
use_samps = paste('S',19:72, sep='')
use_env = c('S','Light_mean','Light_high','Temp_max','Vpd_mean','Vpd_daysatfreq')
Xdata = samples[use_samps,use_env]

# Define tree that each sample comes from
tree = factor(samples[use_samps,'TreeID'])
names(tree) = use_samps

# Extract scores from samples and species only
axes1 = scores(rda0_1, 'sites', choices=1:length(eigenvals(rda0_1)))
axes1_samp = axes1[rownames(Xdata),]
axes1_sp = axes1[sp_samples,]
axes2 = scores(rda0_2, 'sites', choices=1:length(eigenvals(rda0_1)))
axes2_samp = axes2[rownames(Xdata),]
axes2_sp = axes2[sp_samples,]

# Fit vector of species richness
ev1 = envfit(axes1_samp, Xdata, choices=1:4)
ev2 = envfit(axes2_samp, Xdata, choices=1:4)

# Plot ordinations
use_col = rainbow(length(levels(tree)))
names(use_col) = levels(tree)
#color_fact = 1*(rownames(mat1_hel) %in% sp_samples) + 1


pdf(paste(fig_dir, run_name, '_PCA_colbytree_', Tmethod, '.pdf', sep=''), height=9, width=18)
par(mar=c(4,4,3,1))
par(mfrow = c(1,2))

for(use_axes in list(1:2, 3:4)){

# DYE 1
# Define ordination and axes to plot
use_ord = rda0_1

# Calculate percent of variation that each axis explains
ord_sum = summary(use_ord)$cont$importance
pcts = paste(format(ord_sum[2,use_axes]*100, digits=2), '%', sep='')

# Plot
op = ordiplot(use_ord, use_axes, type='n', xlab=paste('PC',use_axes[1],' (',pcts[1], ')', sep=''), 
	ylab=paste('PC',use_axes[2],' (',pcts[2], ')', sep=''), main=names(good_dye)[1], las=1)

# Color samples by trees
for(g in levels(tree)){
	ordihull(op, groups=tree[rownames(mat1_hel)], draw='polygon', col=use_col[g], show.groups=g, alpha=50)
	ordihull(op, groups=tree[rownames(mat1_hel)], draw='line',lwd=2, col=use_col[g], show.groups=g)
}
points(axes1_samp[use_samps, use_axes], pch=21, bg=use_col[tree])

# Add samples from single species
text(axes1_sp[,use_axes], labels=rownames(axes1_sp), cex=.7)

# Add env variables
plot(ev1, choices=use_axes, col='black', add=T, cex=0.9)

# DYE 2
# Define ordination and axes to plot
use_ord = rda0_2

# Calculate percent of variation that each axis explains
ord_sum = summary(use_ord)$cont$importance
pcts = paste(format(ord_sum[2,use_axes]*100, digits=2), '%', sep='')

# Plot
op = ordiplot(use_ord, use_axes, type='n', xlab=paste('PC',use_axes[1],' (',pcts[1], ')', sep=''), 
	ylab=paste('PC',use_axes[2],' (',pcts[2], ')', sep=''), main=names(good_dye)[2], las=1)

# Color samples by trees
for(g in levels(tree)){
	ordihull(op, groups=tree[rownames(mat2_hel)], draw='polygon', col=use_col[g], show.groups=g, alpha=50)
	ordihull(op, groups=tree[rownames(mat2_hel)], draw='line',lwd=2, col=use_col[g], show.groups=g)
}
points(axes2_samp[use_samps, use_axes], pch=21, bg=use_col[tree])

# Add samples from single species
text(axes2_sp[,use_axes], labels=rownames(axes2_sp), cex=.7)

# Add env variables
plot(ev2, choices=use_axes, col='black', add=T, cex=0.9)

}
dev.off()


### Calculate distances between sample and species profiles in ordination space
# Only works for mycobionts
pcdist1 = as.matrix(dist(axes1[,2:15]))[exp_samples, sp_samples]
pcdist2 = as.matrix(dist(axes2[,2:15]))[exp_samples, sp_samples]

pdf(paste(fig_dir, run_name, '_PCdist_samp-sp_',Tmethod,'.pdf', sep=''), height=8, width=20)
layout(matrix(1:4, nrow=1, byrow=T), widths=c(0.3,.3,.3,.1))
par(mar=c(5,5,2,1))

# Presence-Absence community matrix
image(1:Nsp, 1:Nsamp, t(sampXsp_sub), axes=F, xlab='', ylab='', col=c('white','black'), main='Presence/Absence')
axis(1, at=1:Nsp, labels=sp_samples, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Composition based on profiles from Dye1
image(1:Nsp, 1:Nsamp, t(pcdist1/max(pcdist1)), axes=F, xlab='', ylab='', main=names(good_dye)[1],
	col=colorRampPalette(c('black','white'))(10), breaks=seq(0,1,.1))
axis(1, at=1:Nsp, labels=sp_samples, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# Composition based on profiles from Dye2
image(1:Nsp, 1:Nsamp, t(pcdist2/max(pcdist2)), axes=F, xlab='', ylab='', main=names(good_dye)[2],
	col=colorRampPalette(c('black','white'))(10), breaks=seq(0,1,.1))
axis(1, at=1:Nsp, labels=sp_samples, las=2)
axis(2, at=1:Nsamp, labels=exp_samples, las=2)

# BW legend
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
plotColorRamp(c('black','white'), 10, c(.1,.2,.3,.9), labels=seq(0,1,.1))

dev.off()




