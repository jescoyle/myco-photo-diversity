## This script is used to clean up TRFLP profile data which has already been aligned to a size standard.
## It uses a peak size data table exported from PeakScanner2 to:
## 1. Remove peaks from un-used dyes and from the size standard
## 2. Remove undersized peaks due to noise (thresholding)
## 3. Remove peaks from procedural controls


options(stringsAsFactors=F)

# Set data, working, and code directories
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/'
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/TRFLP/photo_TRFLP/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/GitHub//'

setwd(working_dir)

# Read in functions for the analysis of TRFLP data
source('TRFLP_functions.R')

