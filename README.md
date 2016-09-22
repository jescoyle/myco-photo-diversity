# myco-photo-diversity
Contains scripts for the analysis of lichen mycobiont and photobiont communities using T-RFLP profile data. T-RFLP profiles are based on RE digests of fungal and algal specific amplications of portions of ITS. Scripts allow for analysis of multiplexed profiles by different dyes.

TRFLP_functions.R
  Contains functions that are utilized across multiple scripts. Implements fixed percentage, variable percentage, and statistical thresholding of profiles for noise peaks. Includes functions for filtering peaks based on procedural control samples and plotting profiles.

clean_profiles.R
  These scripts take a table of peak sizes exported by PeakScanner and threshold and align peaks to produce profiles for further analysis.
  
analyze_photo_trflp.R
  Compares experimental profiles to one another as well as to expected fragment lengths from in silico RE digests.Conducts ordination of profiles to explore environmental causes of genetic diversity variation.
  

