## This script assess community composition in the eleven different photobiont community data matrices

options(stringsAsFactors=F)

library(reshape2)
library(cluster)
library(plyr)
library(sp)

# Set data, working, and code directories
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Data/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/GitHub/myco-photo-diversity/'
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Figures/'
derived_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Derived_Data/'

setwd(working_dir)

# Load functions
source(file.path(git_dir,'analysis_functions.R', sep=''))

# Set working colors
mycolor = read.csv('../blue2red_10colramp.txt')[9:2,]
mycolor = rgb(mycolor, maxColorValue=255)

# Read in table of info on runs
runs = read.csv('./Analysis/run_names.csv', row.names=1)

# Find community data files
comm_dir = file.path(derived_dir, 'Community_Data_Matrices')
data_files = list.files(comm_dir)

# Read in list of standardized community data
comm_std = sapply(data_files, function(f){
	# Read in data
	this_comm = read.csv(file.path(comm_dir, f), row.names=1, check.names=F)

	# Standardize by total peak area
	this_comm / rowSums(this_comm)
})
names(comm_std) = gsub('.csv', '', names(comm_std))

# Convert NA to zero (NAs are from samples with no peaks)
for(i in 1:length(comm_std)){
	metacomm = comm_std[[i]]
	comm_std[[i]][is.na(metacomm)] = 0
}

# Separate into list for species profiles and sample profiles
sp_std = sapply(comm_std, function(x){
	sp_names = grep('[a-z]+-[a-z]+', rownames(x))
	x[sp_names,]
})

samp_std = sapply(comm_std, function(x){
	S_names = grep('^S[0-9]+$', rownames(x))
	x[S_names,]
})

# Read in environmental data
samples = read.csv(file.path(data_dir, 'samples.csv'))
env = read.csv(file.path(data_dir, 'loggerdata.csv'))
samples = merge(samples, env)
keep_cols = c('SampID','TreeID','Year','BranchPos','Pair','Diam','Angle','Aspect',
	'BryoAbun','Height', 'PosR','PosTheta','Temp_mean','Temp_max','Light_mean','Light_high','Vpd_mean','Vpd_daysatfreq')
samples = samples[,keep_cols]

# Add unique identifier to each branch
samples$BranchID = with(samples, paste(TreeID, BranchPos, sep='.'))

# Add lat/lon to samples
trees = read.csv(file.path(data_dir, 'trees.csv'))
samples = merge(samples, trees[,c('TreeID','Lat','Lon')])

rownames(samples) = paste0('S',samples$SampID)

# Calculate spatial distance matrix between samples
coordinates(samples) = c('Lon','Lat')
proj4string(samples) = CRS("+proj=longlat +units=km")
samples = spTransform(samples, CRS('+proj=utm +zone=17'))
samples = data.frame(samples)
names(samples)[names(samples)=='Lon'] = 'X'
names(samples)[names(samples)=='Lat'] = 'Y'
samples$PosThetaMath = (90 - samples$PosTheta)*(pi/180)
samples$X = samples$X + samples$PosR*cos(samples$PosThetaMath)
samples$Y = samples$Y + samples$PosR*sin(samples$PosThetaMath)
samples$Z = samples$Height # Could add elevation is needed
dist_mat = dist(samples[,c('X','Y','Z')])
dist_mat[dist_mat==0] = 0.5

# Read in fungal community data
lichen_taxa = read.csv(paste0(data_dir, 'taxa.csv'))
rownames(lichen_taxa) = lichen_taxa$TaxonID
myco_comm = read.csv(paste0(data_dir, 'community.csv'))

taxa = unique(myco_comm$TaxonID) # Find all unique taxon IDs
taxa = taxa[!(taxa %in% c('','Fol'))] # Remove taxa that could not be uniquely identified
taxa = taxa[lichen_taxa[taxa,'Lichenized']==1] # Remove all non-lichenized taxa
taxa = taxa[order(taxa)] # Order alphabetically

myco_comm = subset(myco_comm, TaxonID %in% taxa)

# Make a list of species present in whole samples (comm) including and not including species with unknown genus
# This resolves (whenever possible) thalli IDed only to genus
sp_list = lapply(19:72, function(i){
	all_present = unique(subset(myco_comm, SampID==i)$TaxonID)

	these_taxa = lichen_taxa[all_present,]
	unknown = subset(these_taxa, TaxonConcept=='unknown')
	species = subset(these_taxa, TaxonConcept=='species')
	genera = subset(these_taxa, TaxonConcept=='genus')

	sp_genera = unique(species$Genus)
	unique_genera = subset(genera, !(Genus %in% sp_genera))

	new_species = c(species$TaxonID, unique_genera$TaxonID)
	new_species_unk = c(new_species, unknown$TaxonID)

	list(sp = new_species, sp_unk = new_species_unk)
})

# Make mycobiont community data matrix that includes taxa with unknown genus
sampXmyco_unk = sapply(sp_list, function(x){
	sp = x$sp_unk
	taxa %in% sp
})
rownames(sampXmyco_unk) = taxa
colnames(sampXmyco_unk) = paste0('S', 19:72)
sampXmyco_unk = t(sampXmyco_unk)*1
sampXmyco_unk = sampXmyco_unk[,colSums(sampXmyco_unk)>0]

sampXmyco = sapply(sp_list, function(x){
	sp = x$sp
	taxa %in% sp
})
rownames(sampXmyco) = taxa
colnames(sampXmyco) = paste0('S', 19:72)
sampXmyco = t(sampXmyco)*1
sampXmyco = sampXmyco[,colSums(sampXmyco)>0]

# Save
#write.csv(sampXmyco_unk, paste0(data_dir, 'sampXmyco_unk.csv'), row.names=T)
#write.csv(sampXmyco, paste0(data_dir, 'sampXmyco.csv'), row.names=T)

sampXmyco_unk = read.csv(paste0(data_dir, 'sampXmyco_unk.csv'), row.names=1)
sampXmyco = read.csv(paste0(data_dir, 'sampXmyco.csv'), row.names=1)

# Names of photobiont community matrices to analyze
use_mat = paste0('M', 1:10)

# Names of environmental variables to analyze
use_env = c('Light_mean','Light_high','Vpd_mean','Vpd_daysatfreq','Temp_mean')

################################################################################
## How does variation in diversity and community composition vary across space 
## and correspond to variation in the environment?

library(vegan)
library(mgcv)

### Algae ###

## Community composition variation partitioned across spatial hierarchy

ords_list = sapply(use_mat, function(m){

	# Define community matrix
	use_comm = samp_std[[m]]
	use_comm[is.na(use_comm)] = 0
	use_comm = use_comm[, colSums(use_comm)>0]

	# Use Hellinger transformation to standardize
	comm_std = decostand(use_comm, 'hellinger')
	
	# Define sample data
	samp_data = samples[rownames(use_comm),]

	# Define environmental model
	env_std = scale(samp_data[,use_env])

	# Ordination using RDA 
	ord_branch = rda(comm_std ~ samp_data[, 'BranchID'] )
	ord_tree = rda(comm_std ~ factor(samp_data[, 'TreeID']))
	ord_env = rda(comm_std ~ env_std)
	ord_envonly = rda(comm_std ~ env_std + Condition(factor(samp_data[,'TreeID'])))

	list(branch = ord_branch, tree = ord_tree, env = ord_env, env_not_tree = ord_envonly)
})

# R2 in a table
scale_tab = sapply(use_mat, function(m){
	ords = ords_list[,m]
	sapply(c('branch','tree'), function(loc){
		ord = ords[[loc]]
		aov = anova(ord)

		r2 = RsquareAdj(ord)
		c(aov[1,], TotDf=sum(aov$Df), TotVar=sum(aov$Variance), R2adj=RsquareAdj(ord)$adj.r.squared)
	})
}, simplify='array')

write.csv(t(as.data.frame(scale_tab)), file.path(fig_dir, 'Manuscript Results', 'algal_comm_rda_tree-branch_anova.csv'))

pdf(file.path(fig_dir, 'Manuscript Results', 'algal_comm_rda_tree-branch_variance.pdf'), height=4, width=6)
par(mar=c(2,5,1,1))
bp=barplot(as.numeric(scale_tab['R2adj','branch',]), las=1, ylim=c(0,1), border=NA, col=1)
barplot(as.numeric(scale_tab['R2adj','tree',]), las=1, ylim=c(0,1), add=T, col='grey60', border=NA)
mtext('Proportion Variance', 2, 3)
legend('right', c('Tree','Branch'), fill=c('grey60','black'), border=NA, bty='n')
axis(1, at=bp, labels=use_mat, tick=F, line=-1)
abline(v=c(2.5, 4.9), col='grey50', lty=2)
text(1.3, 0.95, runs['M1','RE']) 
text(3.7, 0.95, runs['M3','RE']) 
text(8.5, 0.95, runs['M5','RE'])
dev.off()


## Community composition partitioned between spatial and environmental distance matrices

vp_envtree =  sapply(use_mat, function(m){

	# Define community matrix
	use_comm = samp_std[[m]]
	use_comm[is.na(use_comm)] = 0
	use_comm = use_comm[, colSums(use_comm)>0]

	# Use Hellinger transformation to standardize
	comm_std = decostand(use_comm, 'hellinger')
	
	# Define sample data
	samp_data = samples[rownames(use_comm),]

	# Define environmental model
	env_std = scale(samp_data[,use_env])

	# Partition variance using RDA 
	vp = varpart(comm_std, ~env_std, ~factor(samp_data[,'TreeID']))

	vp
})

# R2 in a table
env_tab = sapply(use_mat, function(m){
	ords = ords_list[,m]

	aov_env = anova(ords[['env']])
	aov_envonly = anova(ords[['env_not_tree']])
	R2adj_env = RsquareAdj(ords[['env']])$adj.r.squared
	R2adj_envonly = RsquareAdj(ords[['env_not_tree']])$adj.r.squared

	stats = c(aov_env[1,], R2adj = R2adj_env, aov_envonly[1,],
		vp_envtree['part',m][[1]]$indfract[,'Adj.R.squared'])
	names(stats)[1:5] = paste(names(stats)[1:5], 'env', sep='_')
	names(stats)[6:9] = paste(names(stats)[6:9], 'envonly', sep='_')
	names(stats)[10:13] = paste('R2adj', c('envonly','both','treeonly','unexplained'), sep='_')
	stats
}, simplify='array')

write.csv(t(env_tab), file.path(fig_dir, 'Manuscript Results', 'algal_comm_rda_env-tree_anova.csv'))

pdf(file.path(fig_dir, 'Manuscript Results', 'algal_comm_rda_env-tree_variance.pdf'), height=4, width=6)
par(mar=c(2,5,1,1))
bp=barplot(1-as.numeric(env_tab['R2adj_unexplained',]), las=1, ylim=c(0,1), border=NA, col=1)
barplot(as.numeric(scale_tab['R2adj','tree',]), add=T, col='grey60', border=NA, axes=F)
barplot(as.numeric(env_tab['R2adj_treeonly',]), add=T, col='grey80', border=NA, axes=F)

mtext('Proportion Variance', 2, 3)
legend('right', c('Environment','Both','Tree'), fill=c('black','grey60','grey80'), border=NA, bty='n')
axis(1, at=bp, labels=use_mat, tick=F, line=-1)
abline(v=c(2.5, 4.9), col='grey50', lty=2)
text(1.3, 0.95, runs['M1','RE']) 
text(3.7, 0.95, runs['M3','RE']) 
text(8.5, 0.95, runs['M5','RE'])
dev.off()


	
# Mantel's correlogram within trees (within 8 m)
mtcor_list = sapply(use_mat, function(m){
	# Define community matrix
	use_comm = samp_std[[m]]
	use_comm[is.na(use_comm)] = 0
	use_comm = use_comm[, colSums(use_comm)>0]

	# Use Hellinger distance
	comm_dist = vegdist(decostand(use_comm, 'hellinger'), 'euclidean')

	# Spatial Distances
	use_dist = as.dist(as.matrix(dist_mat)[rownames(comm_std), rownames(comm_std)])

	# Spatial partioning
	mantel.correlog(comm_dist, use_dist, break.pts=c(0,1,2,3,4,seq(7,551, length.out=7)), cutoff=F)
})

# Plot correlograms
use_pch=c(1,16)
use_col = rainbow(length(use_mat))
names(use_col) = use_mat

pdf(file.path(fig_dir, 'Manuscript Results', 'mantel_correlogram_within_trees.pdf'), height=6, width=6)
par(mar=c(4,4,1,1))
make_plot(xlim=c(0,7), ylim=c(-.1,.4), xlab='Distance (m)', ylab='Correlation (r)')
abline(h=0)
for(m in use_mat){
	dat = mtcor_list[,m]$mantel.res
	sig = dat[1:5,'Pr(corrected)']<0.05
	points(dat[1:5,'class.index'], dat[1:5,'Mantel.cor'], type='o', pch=use_pch[1+as.numeric(sig)], col=use_col[m])

}
dev.off()

pdf(file.path(fig_dir, 'Manuscript Results', 'mantel_correlogram_between_trees.pdf'), height=6, width=6)
par(mar=c(4,4,1,1))
make_plot(xlim=c(0, 551), ylim=c(-.4,.4), xlab='Distance (m)', ylab='Correlation (r)')
abline(h=0)
abline(v=7, lty=2)
for(m in use_mat){
	dat = mtcor_list[,m]$mantel.res
	sig = dat[,'Pr(corrected)']<0.05
	points(dat[,'class.index'], dat[,'Mantel.cor'], type='o', pch=use_pch[1+as.numeric(sig)], col=use_col[m])
}
legend('topright', use_mat, col=use_col, pch=16, ncol=2, bty='n')
dev.off()


## Environmental models of individual strains occuring in at least 10% of samples
	
	

	# Define community matrix
	use_comm = samp_std[[m]]>0
	use_comm[is.na(use_comm)] = 0
	use_comm = use_comm[, colSums(use_comm)>0]

	# Calculate frequency of observation
	use_sp = names(which(colSums(use_comm>0)>(0.1*nrow(use_comm))))
	use_comm = use_comm[,use_sp]
	colnames(use_comm) = paste0('s', colnames(use_comm))
	use_sp = paste0('s', use_sp)

	# Define model data
	mod_data = cbind(samples[rownames(use_comm),use_env], use_comm)

	
	# WORKING HERE- NEED TO DECIDE BEST SPLINE MODEL
	# For each species, model its abundance against each environmental variable
	for(sp in use_sp){
	for(x in use_env){
		mod3 = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=3)), data=mod_data, family=binomial(link='logit'))
		mod6 = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=6)), data=mod_data, family=binomial(link='logit'))



	make_plot(range(mod_data[,x]), c(0,1))
	points(mod_data[,x], mod_data[,sp], col='#0000ff90', pch='|')
	xvals = data.frame(seq(min(mod_data[,x]), max(mod_data[,x]), length.out=100))
	names(xvals)=x
	for(mod in list(mod3, mod6))
	mod_pred = predict(mod, newdata=xvals, type='response', se.fit=T)
	polygon(c(xvals[,1], rev(xvals[,1])), c(mod_pred$fit-mod_pred$se.fit, rev(mod_pred$fit+mod_pred$se.fit)),
		col='#0000ff50', border=NA)
	use_lty = ifelse(summary(mod)$s.table[1,4] < 0.05, 1, 3)
	lines(xvals[,1], mod_pred$fit, col='blue', lwd=2, lty=use_lty)



## Environmental models of RFLP richness


### Mycobiont ###

## Community composition variation partitioned across spatial hierarchy

## Community composition partitioned between spatial and environmental distance matrices


## Environmental models of individual species

## Environmental models of species richness




###############################################################################
### How do diversity and composition of algae and mycobiont relate to one another?

## Richness of mycobionts vs RFLP richness



## Ordination of mycobiont communities with visualization of most abundant algal strains






###############################################################################
### Are mycobiont and algal diversity patterns indicative of environmental constraints?


## Comparison with CAMM



################################################################################
### How is ecological inference influenced by methodology?


## Are the same strains isolated from the mixture of two samples as from each of them separately? (Compare S61 + S62)


## How does PCR affect abundance vs presence of detected strains? (Compare two December runs)


## Do RFLPs from different restrictions enzymes yield different inference? (Compare based on their predicted ability to differentiate based on initial sequences)


## How does inference vary with the method used to remove noise from TRFLPs?


## How does inference vary when more RFLs are binned by using wider base-pair bins to assign strain identities?


###############################################################################
### How frequently are the photobiont strains identified by sequencing also 
### recovered in TRFLP profiles?

## In the community where the DNA was extracted?


## In any community containing the mycobiont species from which the seqeunce was obtained?




#################################################################################
### Code from other files

library(corrplot)
library(vegan)

#### Examine mycobiont communities ####

# How many species per sample?
myco_unk_rich = rowSums(sampXmyco_unk)
myco_unk_rich[order(myco_unk_rich)] # max = 24

myco_rich = rowSums(sampXmyco)
myco_rich[order(myco_rich)] # max = 21 

colSums(sampXmyco)

# How many species overall (not including unID genera)?
length(grep('_', colnames(sampXmyco))) # 57 <- use this value in CAMM
length(grep('_', colnames(sampXmyco_unk))) #70

#### Evaluate which peaks likely correspond to photobiont strains ####

## Which peaks are most abundant in species and community samples?
use_mat = paste0('M', 1:10)

## Calculate total abundance across species and samples
compare_abun = sapply(use_mat, function(m){
	samp_abun = colSums(samp_std[[m]])/nrow(samp_std[[m]])
	sp_abun = colSums(sp_std[[m]])/nrow(sp_std[[m]])
	rbind(samp=samp_abun, sp=sp_abun)
})

# Which peaks encompass 90% of the total abundance from community samples?
top90_abun = sapply(use_mat, function(m){
	tot_abun = compare_abun[[m]]['samp',]
	tot_abun = tot_abun[order(tot_abun, decreasing=T)]
	cum_sum = sapply(1:length(tot_abun), function(i) sum(tot_abun[1:i]))
	cutoff = rev(which(cum_sum - .9 <= 0))[1] + 1
	top90 = tot_abun[1:cutoff]
	
	top90	
})

# Plot
pdf(paste0(fig_dir, 'photo_strain_tot_abun_across_methods_compare_sp.pdf'), height=7, width=10)
par(mar=c(4,4,3,1))
for(m in use_mat){
	if(m=='M1') par(mfrow=c(1,2))
	if(m=='M5') par(mfrow=c(1,3))

	use_data = compare_abun[[m]]

	# Define x-axis limits	
	use_xlim=(max(use_data)+.05)*c(-1,1)
	endpts = trunc(use_xlim*10)*.1
	axis_labs = seq(endpts[1], endpts[2], .1)
	
	# Plot
	plot(0,0, type='n', 
		ylim=c(0,700), xlim=use_xlim,
		ylab='Fragment Length (bp)', xlab='Relative Abundance', axes=F,
		main = paste0(m,': ', paste(runs[m,], collapse=' '))
	)
	abline(h=seq(0,700, 20), col='grey80')
	axis(1, at = axis_labs, labels=abs(axis_labs))
	axis(2, las=1)
	segments(0, as.numeric(colnames(use_data)), use_data['samp',], as.numeric(colnames(use_data)), lend=1)
	segments(0, as.numeric(colnames(use_data)), -use_data['sp',], as.numeric(colnames(use_data)), lend=1, col=2)

	# Add text
	midpt = -endpts[1]/2
	text(use_xlim[1], 710, 'Species Profiles', col=2, pos=4)
	text(use_xlim[2], 710, 'Sample Profiles', col=1, pos=2)
	box()

}
dev.off()

## Plot rank abundance distribution of photobionts across community samples
pdf(paste0(fig_dir, 'photo_strain_tot_abun_across_methods.pdf'), height=7, width=8)
for(m in use_mat){
	use_data = samp_std[[m]]
	use_data[is.na(use_data)] = 0
	samp_abuns = colSums(use_data)/nrow(use_data)
	abun_order = order(samp_abuns, decreasing=T)
	
	tot = samp_abuns[abun_order][1]
	i=1
	while(tot < 0.9){
		i = i + 1
		tot = tot + samp_abuns[abun_order][i]
	}
		
	par(mfrow=c(2,1))
	par(mar=c(4,4,3,1))
	plot(as.numeric(names(samp_abuns)), samp_abuns, type='h', xlim=c(0,700), ylim=c(0, max(samp_abuns)+.05), las=1,
		xlab='Fragment Length (bp)', ylab='Avg. Strain Relative Abun.',
		main = paste0(m,': ', paste(runs[m,], collapse=' ')))
	abline(h = samp_abuns[abun_order][i], col=2)	

	ylab_pts = diff(samp_abuns[abun_order][c(1, i)])*seq(0,1,length.out=i) + samp_abuns[abun_order][1]
	par(mar=c(4,4,1,1))	
	plot(1:length(samp_abuns), samp_abuns[abun_order], type='h', ylim=c(0, max(samp_abuns)+.05),
		axes=F, las=1, xlab='Rank', ylab='Avg. Strain Relative Abun.')
	text(i+1, ylab_pts, names(samp_abuns[abun_order])[1:i], adj=c(-.10,-.10))
	arrows(i+1, ylab_pts, 1:i, samp_abuns[abun_order][1:i], length=.07)

}
dev.off()

## How many photoniont strains are there?

# In each community data matrix?
sapply(compare_abun, function(x) sum(x['samp',]>0)) # in samples
sapply(compare_abun, function(x) sum(x['sp',]>0)) # in thalli

# In the upper 90% abundance of each community data matrix?
sapply(top90_abun, function(x) sum(x>0))

# In the upper 95% abundance
top95_abun = sapply(use_mat, function(m){
	tot_abun = compare_abun[[m]]['samp',]
	tot_abun = tot_abun[order(tot_abun, decreasing=T)]
	cum_sum = sapply(1:length(tot_abun), function(i) sum(tot_abun[1:i]))
	cutoff = rev(which(cum_sum - .95 <= 0))[1] + 1
	tot_abun[1:cutoff]
})
sapply(top95_abun, function(x) sum(x>0))

## How many samples does each strain occur in?
strain_occ = sapply(use_mat, function(m){
	occ = colSums(samp_std[[m]]>0)
	occ = occ[occ>0]
	occ[order(occ, decreasing=T)]
})
sapply(strain_occ, function(x) sum(x>1))

## Define number of photobiont species to be used for CAMM parameters for each community data matrix
CAMM_S = data.frame(all=sapply(compare_abun, function(x) sum(x['samp',]>0)))
CAMM_S$abun90 = sapply(top90_abun, function(x) sum(x>0))
CAMM_S$abun95 = sapply(top95_abun, function(x) sum(x>0))
CAMM_S$occ2 = sapply(strain_occ, function(x) sum(x>1))
write.csv(CAMM_S, file.path(derived_dir, 'compare_strain_richness.csv'), row.names=T)

richvals = unique(as.numeric(as.matrix(CAMM_S)))
richvals[order(richvals)]

## Which photobiont strains are most strongly associated with community variation?

# Unconstrained ordination of photobiont communities
samp_ords = sapply(use_mat, function(m){

	# Define community matrix
	use_comm = samp_std[[m]]
	use_comm[is.na(use_comm)] = 0
	use_comm_bin = use_comm > 0 # Presence/absence

	# Use Hellinger transformation to standardize
	comm_std = decostand(use_comm, 'hellinger')
	comm_std_bin = decostand(use_comm_bin, 'hellinger')

	# Ordination using RDA 
	ord = rda(comm_std)
	ord_bin = rda(comm_std_bin)
	
	list(abun = ord, pres = ord_bin)
})

# Examine each ordination
layout(matrix(1:2*length(use_mat), byrow=T))
for(m in use_mat){
for(type in c('abun','pres')){
	ord = samp_ords[type,m][[1]]
	top90 = names(top90_abun[[m]])
	eigs = eigenvals(ord)
	
	# Variance explained by axes
	var_exp = eigs/sum(eigs)
	cum_var = sapply(1:length(var_exp), function(i) sum(var_exp[1:i]))
	cum_axes = 1:(rev(which(cum_var<.9))[1] + 1)
	
	# Compare to broken-stick model
	# Sometime last axis is greater than broken stick, so we exclude that here
	bstick_axes = 1:(which(bstick(ord) > eigenvals(ord))[1]-1)
	
	# Axes correlated with most abundant strains
	sp_scores = scores(ord, choices=1:length(eigs), display='species')
	top90_axes = apply(sp_scores[top90,], 1, function(x) which(abs(x)==max(abs(x))))
	
	sp = screeplot(ord, bstick=T, npcs=max(cum_axes), las=2, main=paste(m,type))
	abline(v=sp$x[max(bstick_axes)])

	plot_axes = max(bstick_axes) + max(bstick_axes)%%2 # Need an even number of axes
	spnames = rownames(sp_scores)
	sitenames = rownames(samp_std[[m]])

	site_goodness = goodness(ord, display='sites', choices=1:plot_axes, statistic='explained')
	use_col = colorRampPalette(c('white','black'))(10)

	pdf(paste0(fig_dir, 'Photobiont Community RDA/photo_ord_',type,'_',m,'.pdf'), height=4*ceiling(plot_axes/4), width=8)
	par(mfrow=c(ceiling(plot_axes/4),2))
	par(mar=c(4.5,4.5,1,1))
	for(i in seq(1, plot_axes, 2)){
		colorby = cut(site_goodness[,i+1], seq(0,1,.1), include.lowest=T)
		ordiplot(ord, choices=c(i,i+1), type='none', las=1)
		ordihull(ord, choices=c(i,i+1), groups=samples[sitenames, 'TreeID'], display='sites', 
			draw='polygon',label=F, col=1, alpha=50, border='transparent')
		points(ord, choices=c(i,i+1), display='sites', cex=.9, pch=21, col='black', bg=use_col[colorby])
		text(ord, choices=c(i,i+1), 'species', select=spnames %in% top90, cex=1.2)

		these_sp = top90[top90 %in% names(top90_axes[top90_axes %in% c(i, i+1)])]
		if(length(these_sp) >0) text(ord, choices=c(i,i+1), 'species', col=2, cex=1.2, font=2, select=spnames %in% these_sp)	
	}
	dev.off()
	
	# Find species most strongly correlated with each axis
	PC_sp = apply(sp_scores[,bstick_axes], 2, function(x) rownames(sp_scores)[which(x==max(x))])
	important_sp = unique(c(PC_sp, top90))

	# Proportion fit for species
	sp_goodness=goodness(ord, display='species', choices=bstick_axes, statistic='explained')[important_sp,]
	sp_good_byaxis = apply(sp_goodness, 1, function(x) x-c(0,x[1:(length(x)-1)]))

	# Write out table that can be used to fill in strain identities
	df = data.frame(RE = runs[m,'RE'], bp = important_sp, tot_abun = compare_abun[[m]]['samp',important_sp], 
		rank_abun=rank(1-compare_abun[[m]]['samp',])[important_sp], t(sp_good_byaxis))
	df$axis = NA
	df[PC_sp,'axis'] = 1:length(PC_sp)
	write.csv(df, paste0(fig_dir, 'Photobiont Community RDA/important_strains_',m,'_',type,'.csv'), row.names=F)	

	ymax = ceiling(top90_abun[[m]][1]*100)/100

	pdf(paste0(fig_dir,'Photobiont Community RDA/photo_ord_',type,'_',m,'_strains_abun.pdf'), height=4, width=6)
	par(mar=c(4,4,1,1))	
	bp = barplot(t(t(sp_good_byaxis[,top90]) * top90_abun[[m]]), las=2, legend=T, 
		ylim=c(0,ymax), ylab=c('Relative Abundance'), args.legend=list(x='topright', bty='n'))
	points(bp, top90_abun[[m]], pch='-', col=2, lwd=3, lend=1)
	dev.off()
}
}

# For each restriction enzyme, match the important species with the expected strain TRF lengths
exp_len = read.csv(paste0(derived_dir, 'expected_photobiont_strain_length_full.csv'), row.names=1)

REs = unique(runs[,'RE'])
keep_cols = c('RE','bp','tot_abun','rank_abun')

# Compile all important species df into one df
ID_df = data.frame()
for(m in use_mat){
	abun_df = read.csv(paste0(fig_dir, 'Photobiont Community RDA/important_strains_',m,'_abun.csv'))[,keep_cols]
	pres_df = read.csv(paste0(fig_dir, 'Photobiont Community RDA/important_strains_',m,'_pres.csv'))[,keep_cols]

	abun_df$ord_method = 'abun'
	pres_df$ord_method = 'pres'
	df = rbind(abun_df, pres_df)
	
	df$matrix = m
	df = cbind(df, runs[m,c('run','thresh_method')])

	ID_df = rbind(ID_df, df)
}

# Add column giving most likely strain
match_frags = unique(ID_df[,c('RE','bp')])
match_frags$strains = sapply(1:nrow(match_frags), function(i){
	target = match_frags[i,'bp']
	use_lens = exp_len[,match_frags[i,'RE']]

	distances = abs(use_lens-target)
	rownames(exp_len[which(distances==min(distances)),])
})
match_frags$strain_bps = sapply(1:nrow(match_frags), function(i){
	target = match_frags[i,'bp']
	use_lens = exp_len[,match_frags[i,'RE']]

	distances = abs(use_lens-target)
	closest = exp_len[which(distances==min(distances)), match_frags[i,'RE']]
	names(closest) = rownames(exp_len[which(distances==min(distances)),])
	closest
})
match_frags$bp_diff = sapply(1:nrow(match_frags), function(i){
	this_bp = match_frags[i,'bp']
	these_strains = match_frags[i,'strain_bps'][[1]]	
	abs(this_bp - these_strains[1])	
})

# Add species names from which photobiont strains were isolated
strain_names = read.csv(paste0(derived_dir, 'photobiont_strain_names.csv'))

match_frags$species = sapply(1:nrow(match_frags), function(i){
	these_strains = match_frags[i,'strains'][[1]]
	subset(strain_names, Strain_photo %in% these_strains)$SpID
})

match_frags$likely_species = ifelse(match_frags$bp_diff <=5, match_frags$species, NA)

# Save
match_frags_df = match_frags
match_frags_df$strains = sapply(match_frags_df$strains, function(x) paste(unlist(x), collapse=', '))
match_frags_df$strain_bps = sapply(match_frags_df$strain_bps, function(x) paste(unlist(x), collapse=', '))
match_frags_df$species = sapply(match_frags_df$species, function(x) paste(unlist(x), collapse=', '))
match_frags_df$likely_species = sapply(match_frags_df$likely_species, function(x) paste(unlist(x), collapse=', '))
write.table(match_frags_df, paste0(fig_dir, 'Photobiont Community RDA/match_RFLs.txt'), sep='\t', row.names=F)



#### WORKING HERE EXPLORING WHICH SPECIES ARE ASSOCIATED WITH THE MOST IMPORTANT PHOTOBIONT STRAINS


## For each community sample generate a list of possible strains from each community matrix
strain_arr = array(NA, dim=c(nrow(sampXmyco), nrow(exp_len), length(comm_std)), 
	dimnames=list(SampID=rownames(sampXmyco), strain=rownames(exp_len), matrix=names(comm_std)))

for(m in names(comm_std)){
	use_lens = exp_len[,runs[m,'RE']]
	names(use_lens) = rownames(exp_len)
	
	this_mat = samp_std[[m]]

	# Drop columns without peaks
	this_mat = this_mat[,colSums(this_mat)>0]
	

	for(p in names(use_lens)){
		these_cols = abs(use_lens[p] - as.numeric(names(this_mat)))<= 5
		if(sum(these_cols)==0){
			strain_arr[,p,m] = 0
		} else {
			if(sum(these_cols) > 1){
				these_abuns = colSums(this_mat[,these_cols])
				this_peak = names(which(these_abuns==max(these_abuns)))
			} else {
				this_peak = names(this_mat)[these_cols]
			}
			strain_arr[,p,m] = this_mat[,this_peak]
		}		
	}
}
strain_molten = melt(strain_arr)


## Across methods, are strains found repeatedly?

Nmat = length(use_mat)
Nsamp = nrow(sampXmyco) 


image(1:Nmat, 1:Nsamp, t(strain_arr[,1,use_mat])>0, col=c('white','black'), axes=F, xlab='', ylab='')
axis(3, at=1:Nmat, labels=use_mat, tick=F, line=-1)
axis(2, at=1:Nsamp, labels=rownames(sampXmyco), tick=F, line=-.8,  las=1)
axis(1, at=c(1.5, 3.5, 7.5), labels=runs[c(1,3,5), 'RE'], line=-1, tick=F, col=2)
abline(v=(0:Nmat)+0.5, col='grey80')
abline(h=(0:Nsamp)+0.5, col='grey80')
abline(v=c(2.5, 4.5), col=1, lwd=2)









### THIS DOESN'T MAKE SENSE B/C ABUNDANT SPECIES POOLS ACROSS SAMPLES
# For each community matrix, compare known mycobiont species composition to likely species in photobiont strains
for(m in use_mat){
for(type in c('abun','pres')){
	
	# Compile list of species from which most abundance photobionts may have come
	imp_sp = subset(ID_df, matrix==m & ord_method==type)
	imp_sp = merge(imp_sp, match_frags)
	splist = na.omit(unique(unlist(imp_sp$likely_species)))

	# Get list of mycobiont t
	


}

# Add info
#ID_df = merge(ID_df, match_frags)







# Check whether these peaks were the ones expected.

## Hierarchical clustering of species profiles

sp_dists = dist(comm_std)
cl = hclust(sp_dists, method='ward.D2')
plot(cl)





#### Variation in communities with environment vs spatial distance####

# Calculate species richness for each photobiont community matrix
photo_rich = sapply(samp_std, function(x) rowSums(x>0))
photo_rich[is.na(photo_rich)] = 0

# Correspondence of trflp methods
cor_mat = cor(photo_rich, method='spearman')
pdf(paste0(fig_dir, 'rho_photobiont_richness_across_methods.pdf'), height=5, width=5)
corrplot.mixed(cor_mat, order='AOE', tl.pos='d', lower='number',upper='ellipse')
dev.off()

# Calculate species richness for mycobiont community matrices
myco_rich = data.frame(species=rowSums(sampXmyco), all=rowSums(sampXmyco_unk))

# Calculate total richness across all samples
photo_rich_tot = sapply(samp_std, function(x) sum(colSums(x, na.rm=T)>0))
myco_rich_tot = c(species=sum(colSums(sampXmyco)>0), all=sum(colSums(sampXmyco_unk)>0))

## Plot richness vs env for each method

use_env = c('Light_mean','Light_high','Vpd_mean')
use_env = c('Light_mean','Light_high','Vpd_mean','Vpd_daysatfreq','Temp_mean','Temp_max')
use_col = rainbow(9)

pdf(paste0(fig_dir, 'compare_richness_vs_env_across_methods.pdf'), height=6, width=10)
layout(matrix(1:6, ncol=3, byrow=F))
par(mar=c(4,4,1,1))

for(xvar in use_env){
	xdata = samples[rownames(photo_rich),xvar]
	
	plot(xdata, rep(0, length(xdata)), type='n', ylim=c(0,45), las=1, xlab=xvar, ylab='Photobiont Richness')

	for(i in 1:9){
		ydata = photo_rich[,paste0('M',i)]
		points(xdata, ydata, pch=16, col=use_col[i])
		abline(lm(ydata ~ xdata), col=use_col[i])
	}
	legend('topright', paste0('M',1:9), pch=16, col=use_col, ncol=2)

	plot(xdata, rep(0, length(xdata)), type='n', ylim=c(0,25), las=1, xlab=xvar, ylab='Mycobiont Richness')
	for(i in 1:2){
		ydata = myco_rich[,i]
		points(xdata, ydata, pch=16, col=i)
		abline(lm(ydata ~ xdata), col=i)
	}
	legend('topright', colnames(myco_rich), pch=16, col=1:2)
}
dev.off()



## Which strains are most associated with overall community variation vs to environmental variation?



## Calculate RDA of community composition
env_data = samples[rownames(samp_std[[1]]), use_env]

photo_rda_pres = sapply(paste0('M', 1:11), function(m){
	sapply(use_env, function(x){
		use_comm = samp_std[[m]]
		use_comm[is.na(use_comm)] = 0
		unlist(calc_rda(use_comm, env_data[,x], binary=T))		
	}, simplify=T)
}, simplify='array')
names(dimnames(photo_rda_pres)) = c('Statistic','Env','Matrix')

photo_rda_abun = sapply(paste0('M', 1:11), function(m){
	sapply(use_env, function(x){
		use_comm = samp_std[[m]]
		use_comm[is.na(use_comm)] = 0
		unlist(calc_rda(use_comm, env_data[,x], binary=F))		
	}, simplify=T)
}, simplify='array')
names(dimnames(photo_rda_abun)) = c('Statistic','Env','Matrix')

myco_rda_pres = sapply(list(species=sampXmyco, all=sampXmyco_unk), function(m){
	sapply(use_env, function(x){
		use_comm = m
		use_comm[is.na(use_comm)] = 0
		unlist(calc_rda(use_comm, env_data[,x], binary=T))		
	}, simplify=T)
}, simplify='array')
names(dimnames(myco_rda_pres)) = c('Statistic','Env','Matrix')

molten = melt(photo_rda_pres)
molten$Matrix = factor(molten$Matrix, levels=paste0('M', 1:11))
photo_rda = dcast(molten, Env + Matrix ~ Statistic)

molten = melt(myco_rda_pres)
myco_rda = dcast(molten, Env + Matrix ~ Statistic)

bp_data = dcast(molten, Matrix ~ Env, subset=.(Statistic=='R2'))
rownames(bp_data) = bp_data[,1]
bp_data = as.matrix(bp_data[,-1])

# only run the figure that pertains to the molten dataframe you are using
pdf(file.path(fig_dir, 'compare_RDA_R2_photo_across_methods.pdf'), height=4, width=10)
bp = barplot(bp_data, beside=T, las=1, ylab=expression(RDA~~R^2), space=c(0.2, 1), ylim=c(0,.1))
text(bp, bp_data, labels=1:11, pos=3)
dev.off()

pdf(file.path(fig_dir, 'compare_RDA_R2_myco_across_methods.pdf'), height=4, width=10)
bp = barplot(bp_data, beside=T, las=1, ylab=expression(RDA~~R^2), space=c(0.2, 1), ylim=c(0,0.05))
text(bp, bp_data, labels=c('species','all'), pos=3)
dev.off()


################################################
### Compare species composition to species photobiont profiles


# Match species peaks to sample peaks to see whether species are in samples or match other species
# ONLY RUN ONCE THE RE-LOAD BELOW
spINcomm = sapply(names(samp_std), function(m){
	print(m)
	metacomm = rbind(sp_std[[m]],samp_std[[m]])
	species = sp_std[[m]]
	
	# Make blank array to hold matches
	matching = array(NA, dim=c(nrow(species), nrow(metacomm),3),
		dimnames=list(SpID=rownames(species), SampID=rownames(metacomm), Statistic=c('Pct_match','N_match','N_peaks')))

	# Calculate matches for each species profile in each sample profile
	for(i in rownames(species)){
	for(j in rownames(metacomm)){
		matching[i,j,] = match_peaks(species[i,], metacomm[j,])
	}}

	matching
})
names(spINcomm) = names(samp_std)

# Save list of arrays
#save(spINcomm, file=paste0(derived_dir,'match_species_profiles_to_comm.RData'))
load(paste0(derived_dir,'match_species_profiles_to_comm.RData'))

# Examine which peaks species have in common
pdf(paste0(fig_dir, 'pct_match_species profiles.pdf'), height=6, width=11)
par(mar=c(4.5,4.5,3,1))
par(mfrow=c(1,2))
for(m in paste0('M', 1:11)){
	spnames = dimnames(spINcomm[[m]])$SpID
	image(spINcomm[[m]][spnames,spnames,'Pct_match'], main=paste(runs[m,], collapse=' '), 
		axes=F, col=colorRampPalette(c('white','black'))(10))
	axis(1, at=seq(0,1, length.out=length(spnames)), labels=spnames, las=2)
	axis(2, at=seq(0,1, length.out=length(spnames)), labels=spnames, las=2)

	sp_dist = dist(t(spINcomm[[m]][spnames,spnames,'Pct_match']))
	sp_clust = hclust(sp_dist, method='average')
	plot(sp_clust, xlab='')
}
dev.off()

## Compare known community composition to species profiles
# Note: 'lec-sp2' was IDed as 'lec-str'.
# Photobionts: have 'lec-str' only so leave name in community matrix alone
# Photobionts: have 'usn-sp2' which should be 'usn-str' in community matrix
# Both were sequenced and had different photobiont strains to don't change names- this means that we won't be searching for photobiont from usn-sp2

# Define metacommunity matrix
sampXsp = sampXmyco_unk

# Convert column names and subset to species for which we have profiles
colnames(sampXsp) = tolower(sub('_', '-', colnames(sampXsp)))

# Find difference between proportion of a species' profile that is in a sample and whether the species is in the sample (0,1)
# Will be 0 when full profile in sample and species in sample
# Will be positive when species in sample and some of profile in sample
# Will be negative when species not in sample but some of profile is in sample
# Magnitude gives dissimilarity (0, 1)
spsamp_dist = sapply(names(spINcomm), function(m){
	#print(m)
	species = dimnames(spINcomm[[m]])$SpID
	species = species[!(species %in% c('lec-sp2','usn-sp2'))]
	sampXsp[,species] - t(spINcomm[[m]][species, rownames(sampXsp),'Pct_match'])
})

# Plot

# Mycobiont presence/absence matrix
Nsp = ncol(sampXsp)
Nsamp = nrow(sampXsp)

pdf(paste0(fig_dir, 'mycobiont_community_matrix.pdf'), height=9.5, width=11.5)
image(1:Nsp, 1:Nsamp, t(sampXsp), axes=F, xlab='', ylab='', col=c('transparent','black'), main='Presence/Absence')
abline(v=(0:Nsp)+0.5, col='grey80')
abline(h=(0:Nsamp)+0.5, col='grey80')
axis(1, at=1:Nsp, labels=colnames(sampXsp), las=2)
axis(2, at=1:Nsamp, labels=rownames(sampXsp), las=2)
dev.off()

use_col = c(mycolor[1:4], '#FFFFFF00', mycolor[5:8])

pdf(paste0(fig_dir, 'compare_sp_samp_profiles.pdf'), height=20, width=12)
layout(matrix(1:6, nrow=3, byrow=T))
par(mar=c(5,5,2,1))

for(m in paste0('M', 1:11)){
	use_dist = spsamp_dist[[m]]

	species = 

	# Define matrix layout based on metacommunity matrix
	Nsp = ncol(use_dist)
	Nsamp = nrow(use_dist)

	# Dissimilarity from expected
	image(1:Nsp, 1:Nsamp, t(use_dist), axes=F, xlab='', ylab='', main=paste0(runs[m,], collapse=' '),
		col=colorRampPalette(use_col)(20), breaks=seq(-1,1,.1))
	abline(v=(0:Nsp)+0.5, col='grey80')
	abline(h=(0:Nsamp)+0.5, col='grey80')
	axis(1, at=1:Nsp, labels=colnames(use_dist), las=2)
	axis(2, at=1:Nsamp, labels=rownames(use_dist), las=2)
}

# Color legend
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
plotColorRamp(use_col, 20, c(.5,.2,.6,.9), labels=seq(-1,1,.2))
dev.off()




## WORKING HERE














