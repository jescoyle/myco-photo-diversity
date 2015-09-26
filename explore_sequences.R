## This script analyzes sequence data for the Mycobiont - Photobiont Diversity Project

options(stringsAsFactors=F)

# Install Bioconductor packages
#source("http://bioconductor.org/biocLite.R")
#nbiocLite('REDseq') # 'Biostrings'


library(seqinr)
library(Biostrings)
browseVignettes(package = "Biostrings")

seq_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Sequences/al1700f-its4'
setwd(seq_dir)

# Find sequence files in the sequences directory
seq_filenames = list.files(seq_dir, pattern='*.TXT', full.name=F)

seqs = sapply(seq_filenames, function(x) read.fasta(x))

taxa = sapply(strsplit(seq_filenames, '-'), function(x) x[1])

seqs_bio = DNAStringSet(sapply(seqs, function(x) toupper(c2s(x))), use.names=F)
names(seqs_bio) = taxa

# Calculate length of each sequence
seqs_len = sapply(seqs, length)
names(seqs_len) = taxa

## Simple pairwise alignment of each sequence to determine which are more similar
# Set substitution penalties matrix
sigma = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = F) # includes ambiguous nucleotides

# Alignment

use_pattern = sapply(seqs, function(x) toupper(c2s(x)))
use_subject = toupper(c2s(seqs[[1]]))

# Make matrix of pairwise alignment scores
use_pattern = sapply(seqs, function(x) toupper(c2s(x)))
score_mat = sapply(seqs, function(x){
	use_subject = toupper(c2s(x))
	pairwiseAlignment(subject=use_subject, pattern=use_pattern, substitutionMatrix=sigma,
	gapOpening=-2, gapExtension=-.5, scoreOnly=T)
})

rownames(score_mat) = taxa
colnames(score_mat) = taxa

# Standardize by max score (on diag)
score_mat_std = score_mat/diag(score_mat)

library(corrplot)
corrplot(score_mat_std, method='square')

# Simple hierarchical clustering
dists = as.dist(1-score_mat_std)
cluster <- hclust(dists, method="ward.D2")
plot(cluster)

## Align groups
set1a = seqs_bio[c('lec_cae','usn_sp2','bue_eru','usn_sp1')]
set1b = seqs_bio[c('pmo_sub','lec_alb','lec_sp1','och_afr','pyr_var','lep_sp1','per_sp3')]
set2a = seqs_bio[c('per_tra','pyx_sub','mye_aur','phy_pum','phy_ste')]
set2b = seqs_bio[c('pmo_per','htr_liv','usn_str','pmo_sp1','pmo_sp2','art_sus','pmo_sui','pmo_ult')]
set2c = seqs_bio[c('can_cro','can_car','mye_obs','can_tex','pmo_ret','pps_min')]
set3a = seqs_bio[c('fla_cap','phy_mil','phy_ame','pun_rud')]
set3ai = seqs_bio[c('het_obs','pmo_hyp','pmo_hyp2')]
set3b = seqs_bio[c('lec_str','unid_cr3')]
setTrent = seqs_bio[c('art_rub','gra_scr')]

writeXStringSet(set1a, 'group_1a.fasta')
writeXStringSet(set1b, 'group_1b.fasta')
writeXStringSet(set2a, 'group_2a.fasta')
writeXStringSet(set2b, 'group_2b.fasta')
writeXStringSet(set2c, 'group_2c.fasta')

writeXStringSet(set3a, 'group_3a.fasta')
writeXStringSet(set3ai, 'group_3ai.fasta')
writeXStringSet(set3b, 'group_3b.fasta')
writeXStringSet(setTrent, 'group_Trent.fasta')

writeXStringSet(seqs_bio, 'jan2015_seqs.fasta')



####### Compare new sequences to established strains

# Read in strains
strains = read.fasta('../photobiont_types/photobiont_sequences.fas')

# Read in new taxon sequence
this_sp = 'per_sp4'

this_seq = read.fasta(paste(this_sp,'pS','al1700f_its4.TXT', sep='-'), 
	as.string=T, seqonly=T)

# Compare sequence
# Set substitution penalties matrix
sigma = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = F) # includes ambiguous nucleotides
use_pattern = sapply(strains, function(x) toupper(c2s(x)))
lineup = pairwiseAlignment(subject=unlist(this_seq), pattern=use_pattern, substitutionMatrix=sigma,
	gapOpening=-2, gapExtension=-.5, scoreOnly=T)
which(lineup==max(lineup))

##############################################################
### Find cut sites and digest sequences

library(Biostrings)
library(FNN) # search nearest neighbors knn.dist

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
setwd(working_dir)

photoseqs = readDNAStringSet('./Sequences/photobiont_types/photobiont_sequences_ITS1-ITS2_nodups.fas')

# Expected fragment lengths (including 19bp ITS1 and 20bp ITS2 primers)
its12_frags = data.frame(strain=names(photoseqs),length=width(photoseqs)+39)
its12_frags[order(its12_frags$length),]
# Could probably distinguish 7-8 strains from length alone

# Read in table of fragment length from REPK (http://rocaplab.ocean.washington.edu/tools/repk)
# These do not include length added by primer 
frag_tab_fw = read.table('./TRFLP/repk_its1-its2.txt', header=T, sep='\t', row.names=1)
frag_tab_rv = read.table('./TRFLP/repk_its2-its1.txt', header=T, sep='\t', row.names=1)

fw_unifrag = apply(frag_tab_fw, 2, function(x) length(unique(x)))
rv_unifrag = apply(frag_tab_rv, 2, function(x) length(unique(x)))

unifrags = data.frame(re = colnames(frag_tab_fw), fw = fw_unifrag, rv = rv_unifrag)

max(unifrags$fw)
max(unifrags$rv)
subset(unifrags, fw>17|rv>17)
subset(unifrags, fw>17)

# Which single enzyme gives maximally dispersed fragment lengths?
re_disp = apply(frag_tab, 2, function(x){
	D = knn.dist(x, k=1)
	#D = rowMeans(D)
	lowD = D[D<=quantile(D,.5)]
	#D = D[D>5]
	mean(lowD)
})
use_re = names(re_disp[re_disp>0])

re_disp[order(re_disp)]

# Which pair of enzymes gives maximally dispersed fragment lengths in two dimensions?
re_dmat = matrix(NA, ncol(frag_tab), ncol(frag_tab))
colnames(re_dmat) = colnames(frag_tab)
rownames(re_dmat) = colnames(frag_tab)
for(i in 1:nrow(re_dmat)){
for(j in i:ncol(re_dmat)){
	coords = frag_tab[,c(i,j)]
	#D = dist(coords, method='maximum')
	D = knn.dist(coords, k=1)
	#D = rowMeans(D)
	lowD = D[D<=quantile(D,.5)]
	re_dmat[i,j] = mean(lowD)
}}

hist(re_dmat[upper.tri(re_dmat)])
which(re_dmat >5.5, arr.ind=T)

re_dmat[c(81,93,139,113,113),c(93,113,215,265,283)]

# Plot fragments length profiles under two REs
# Choose enzymes and side considered 5'
re1 = 'Bst4CI'
dir1 = 'rv'
re2 = 'Eco32I'
dir2 = 'fw'

# Plot a single pair
use_x = frag_tab[,paste(re1,dir1,sep='_')]
use_y = frag_tab[,paste(re2,dir2,sep='_')]

pdf('RE AgsI-rv HpyCH4V-fw.pdf', height=6, width=6)

par(mar=c(4,4,.5,.5))
plot(use_x, use_y, type='n', xlab=paste(re1,dir1), ylab=paste(re2,dir2), las=1,
	xlim=c(0, max(width(photoseqs))),ylim=c(0, max(width(photoseqs))))
text(use_x, use_y, labels=rownames(frag_tab))
points(use_x,rep(0,length(use_x)), pch=3)
points(rep(0,length(use_y)),use_y, pch=3)
points(width(photoseqs),rep(-10,length(photoseqs)), pch=3, cex=.5, col=2)
points(rep(-10,length(photoseqs)),width(photoseqs), pch=3, cex=.5, col=2)
axis(1, at=seq(5,500,5), labels=F, tcl=-.3)
axis(2, at=seq(5,500,5), labels=F, tcl=-.3)

dev.off()

## Plot top two pairs
pairs = list(c('BsiSI_fw','BssKI_fw'), c('Hpy188III_fw','BssKI_fw'))

pdf('RE top 2 pairs.pdf', height=5, width=10)
par(mfrow=c(1,2))
par(mar=c(4,4,.5,.5))
for(x in pairs){
	use_x = frag_tab[,x[1]]
	use_y = frag_tab[,x[2]]

	plot(use_x, use_y, type='n', xlab=x[1], ylab=x[2], las=1,
		xlim=c(0, max(width(photoseqs))),ylim=c(0, max(width(photoseqs))))
	text(use_x, use_y, labels=rownames(frag_tab))
	points(use_x,rep(0,length(use_x)), pch=3)
	points(rep(0,length(use_y)),use_y, pch=3)
	points(width(photoseqs),rep(-10,length(photoseqs)), pch=3, cex=.5, col=2)
	points(rep(-10,length(photoseqs)),width(photoseqs), pch=3, cex=.5, col=2)
	axis(1, at=seq(5,500,5), labels=F, tcl=-.3)
	axis(2, at=seq(5,500,5), labels=F, tcl=-.3)
}
dev.off()

## Make lists of undifferentiable groups under each pair
this_pair = pairs[[2]]
D = as.matrix(dist(frag_tab[,this_pair], method='maximum'))
1*(D<5)


## Make plot of strains under two restriction scenarios

# Put tables together
names(frag_tab_fw) = paste(names(frag_tab_fw), 'fw', sep='_')
names(frag_tab_rv) = paste(names(frag_tab_rv), 'rv', sep='_')
frag_tab = cbind(frag_tab_fw, frag_tab_rv)

# Do for all pairs that can distinguish at least 15 unique lengths
use_fw = paste(subset(unifrags, fw>16)$re, 'fw', sep='_')
use_rv = paste(subset(unifrags, rv>16)$re, 'rv', sep='_')
use_re = c(use_fw, use_rv)

pdf('RE lengths most dipersed.pdf', height=6, width=6)
for(i in 1:length(use_re)){
for(j in i:length(use_re)){
	use_x = frag_tab[,use_re[i]]
	use_y = frag_tab[,use_re[j]]
	
	par(mar=c(4,4,.5,.5))
	plot(use_x, use_y, type='n', xlab=use_re[i], ylab=use_re[j], las=1,
		xlim=c(0, max(width(photoseqs))),ylim=c(0, max(width(photoseqs))))
	#text(use_x, use_y, labels=rownames(frag_tab))
	points(use_x, use_y, pch=4)
	points(use_x,rep(0,length(use_x)), pch=3)
	points(rep(0,length(use_y)),use_y, pch=3)
	points(width(photoseqs),rep(-10,length(photoseqs)), pch=3, cex=.5, col=2)
	points(rep(-10,length(photoseqs)),width(photoseqs), pch=3, cex=.5, col=2)
	axis(1, at=seq(5,500,5), labels=F, tcl=-.3)
	axis(2, at=seq(5,500,5), labels=F, tcl=-.3)
}}
dev.off()



# Save expected length table for pair that I decided to use: Hpy188III BssKI
write.csv(frag_tab_fw[,c('Hpy188III','BssKI')], file='expected_photobiont_strain_lengths.csv', row.names=T)

## Generate profiles for control communities

M1 = paste('p', c(1,5,6,9,12,19,24,30), sep='')
M2 = paste('p', c(5,6,12,19), sep='')
M3 = paste('p', c(1,3,5,12), sep='')
M4 = paste('p', c(1,6,19,34), sep='')

frag_tab_fw[M1,c('Hpy188III','BssKI')] + 39 # for primers


## Mycobiont RE choice

mycoseqs = readDNAStringSet('./Sequences/58m2-nlb4/mycobiont_sequence.fas')

strain_names = paste('m', 1:length(mycoseqs), sep='')
match_names = data.frame(seq_label = names(mycoseqs), strainID = strain_names)
rownames(match_names)=match_names$seq_label

write.csv(match_names, file='./Sequences/mycobiont_strain_names.csv', row.names=F)

names(mycoseqs) = strain_names
writeXStringSet(mycoseqs, file='./Sequences/58m2-nlb4/mycobiont_strains.fas')

frag_tab = read.table('./TRFLP/repk_58m2-nlb4.txt', sep='\t', header=T, row.names=1)


# Pair of enzymes giving maximally dispersed fragments in two dimensions (from code above)
hist(re_dmat[upper.tri(re_dmat)])
which(re_dmat >3.8, arr.ind=T)

# Plot likely pairs

re_dmat[c(11,81,135,137),c(11,81,135,137)]

use_re = colnames(frag_tab)[c(11,81,135,137)]

pdf('RE lengths most dipersed mycobiont.pdf', height=6, width=6)
for(i in 1:length(use_re)){
for(j in i:length(use_re)){
	use_x = frag_tab[,use_re[i]]
	use_y = frag_tab[,use_re[j]]
	
	par(mar=c(4,4,.5,.5))
	plot(use_x, use_y, type='n', xlab=use_re[i], ylab=use_re[j], las=1,
		xlim=c(0, max(width(mycoseqs))),ylim=c(0, max(width(mycoseqs))))
	#text(use_x, use_y, labels=rownames(frag_tab))
	points(use_x, use_y, pch=4)
	points(use_x,rep(0,length(use_x)), pch=3)
	points(rep(0,length(use_y)),use_y, pch=3)
	points(width(mycoseqs),rep(-10,length(mycoseqs)), pch=3, cex=.5, col=2)
	points(rep(-10,length(mycoseqs)),width(mycoseqs), pch=3, cex=.5, col=2)
	axis(1, at=seq(5,600,5), labels=F, tcl=-.3)
	axis(2, at=seq(5,600,5), labels=F, tcl=-.3)
}}
dev.off()


pdf('RE lengths with Hpy188III mycobiont.pdf', height=6, width=6)
for(j in 1:length(use_re)){
	use_x = frag_tab[,'Hpy188III']
	use_y = frag_tab[,use_re[j]]
	
	par(mar=c(4,4,.5,.5))
	plot(use_x, use_y, type='n', xlab='Hpy188III', ylab=use_re[j], las=1,
		xlim=c(0, max(width(mycoseqs))),ylim=c(0, max(width(mycoseqs))))
	#text(use_x, use_y, labels=rownames(frag_tab))
	points(use_x, use_y, pch=4)
	points(use_x,rep(0,length(use_x)), pch=3)
	points(rep(0,length(use_y)),use_y, pch=3)
	points(width(mycoseqs),rep(-10,length(mycoseqs)), pch=3, cex=.5, col=2)
	points(rep(-10,length(mycoseqs)),width(mycoseqs), pch=3, cex=.5, col=2)
	axis(1, at=seq(5,600,5), labels=F, tcl=-.3)
	axis(2, at=seq(5,600,5), labels=F, tcl=-.3)
}
dev.off()


## Plot REs used (Hpy188III and AccII)
pdf('./TRFLP/RE lengths with Hpy188III and BstUI(AccII) mycobiont labeled.pdf', height=6, width=6)
	use_x = frag_tab[,'Hpy188III']
	use_y = frag_tab[,'AccII']
	
	par(mar=c(4,4,.5,.5))
	plot(use_x, use_y, type='n', xlab='Hpy188III', ylab='BstUI', las=1,
		xlim=c(0, max(width(mycoseqs))),ylim=c(0, max(width(mycoseqs))))
	text(use_x, use_y, labels=substring(rownames(frag_tab),2))
	#points(use_x, use_y, pch=4)
	points(use_x,rep(0,length(use_x)), pch=3)
	points(rep(0,length(use_y)),use_y, pch=3)
	points(width(mycoseqs),rep(-10,length(mycoseqs)), pch=3, cex=.5, col=2)
	points(rep(-10,length(mycoseqs)),width(mycoseqs), pch=3, cex=.5, col=2)
	axis(1, at=seq(5,600,5), labels=F, tcl=-.3)
	axis(2, at=seq(5,600,5), labels=F, tcl=-.3)
dev.off()


unifrags = apply(frag_tab, 2, function(x) length(unique(x)))
max(unifrags)
unifrags[unifrags>35]

use_re = names(unifrags[unifrags>35])



## Generate known profiles
profiles = frag_tab[,c('Hpy188III', 'AccII')]
write.csv(profiles,'./Analysis/Derived_Data/expected_mycobiont_strain_lengths.csv', row.names=T)

C1 = paste('m', c(1,5,9,12,15,16,21,30,32,35,37,39,40,47) , sep='')
C2 = paste('m', c(1,5,12,21,30,37,40,47), sep='')

use_sp = C2

use_x = frag_tab[use_sp,'Hpy188III']
use_y = frag_tab[use_sp,'AccII'] # Same as BstUI

pdf('./TRFLP/mycobiont control community2.pdf', height=6, width=6)	
par(mar=c(4,4,.5,.5))
plot(use_x, use_y, type='n', xlab='Hpy188III', ylab='BstUI', las=1,
	xlim=c(0, max(width(mycoseqs))),ylim=c(0, max(width(mycoseqs))))

text(use_x, use_y, use_sp, col=2)
points(use_x,rep(0,length(use_x)), pch=3)
points(rep(0,length(use_y)),use_y, pch=3)

axis(1, at=seq(5,600,5), labels=F, tcl=-.3)
axis(2, at=seq(5,600,5), labels=F, tcl=-.3)
dev.off()



### T-RFLP species

sp = c('Bue_eru','Usn_sp2','Lec_sp1','Pmo_sub','Art_sus','Pmo_ult','Pmo_sui','Pmo_per','Htr_liv',
	'Usn_str','Can_cro','Can_car','Mye_obs','Can_tex','Pmo_ret','Pps_min','Per_tra','Pyx_sub','Mye_aur','Phy_pum',
	'Phy_ste','Fla_cap','Phy_mil','Phy_ame','Pun_rud','Pmo_hyp','Het_obs','Lec_str','Art_rub','Gra_scr',
	'Cnd_con','Mar_pol','Per_sp2','Gya_buc','Ope_vul','Sco_chl','Cdl_eff','Cat_nig','Per_tet')


sp[order(sp)]

## Look at certain species
use_re = c('Hpy188III','AccII')
frag_tab[paste('m',24:30, sep=''),use_re]

match_names


###################################################
### Old Code
# Choose enzymes
enz_fwd = c('BshFI','BsiSI','CviJI', 'SetI') # 
enz_rev = c('AfiI', 'FaiI', 'AgsI', 'SetI') # ,

enz = unique(c(enz_fwd, enz_rev))

frag_fwd = frag_tab[,enz]
frag_rev = frag_tab_rv[,enz]

apply(frag_fwd, 2, function(x) length(unique(x)))
apply(frag_rev, 2, function(x) length(unique(x)))

use_rev = frag_rev[,c('AgsI','AfiI')]
use_rev = use_rev[order(use_rev$AgsI, use_rev$AfiI),]
use_mins = apply(use_rev, 1, min)

# Find minimum frag lengths
lens_fwd = apply(frag_fwd, 1, min)
lens_fwd[order(lens_fwd)]

lens_rev = apply(frag_rev, 1, min)
lens_rev[order(lens_rev)]

# Make table of peak pairs
peakpairs = data.frame(strain=rownames(frag_tab),lens_fwd,lens_rev)
peakpairs[order(peakpairs$lens_fwd, peakpairs$lens_rev),]

# Make table of peak pairs with 2 enzymes
best_fwd = 'SetI'
best_rev = 'AgsI'

fwd_frag = frag_tab[,best_fwd]
rev_frag = frag_tab_rv[,best_rev]

rev_frag[order(rev_frag)]
fwd_frag[order(fwd_frag)]

pairs = data.frame(strain=rownames(frag_tab), fwd_frag, rev_frag)

pair_tab = table(pairs$fwd_frag, pairs$rev_frag)





