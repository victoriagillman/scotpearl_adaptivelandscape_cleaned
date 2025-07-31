### Using plink on VP files- filtering for missing and linkage disequilibrium
#Code modified from Sam Beck :)

# Housekeeping!
rm(list=ls())

# A brief bit of PLINK command line coding
#Set wd to be where .bed .bim .fam files are
setwd("./data/plink_VP/")
getwd()
# Use system coding to talk to plink- check if its working
system("cmd.exe /c dir")
system("./plink.exe")

# First need to make ped and map
system("./plink --bfile PM_Dec4 --allow-extra-chr --recode tab --out pearlmusselfiles")
## make another ped file with recode to get a .raw file
system("./plink --file pearlmusselfiles --allow-extra-chr --recode A --out pearlmusselfiles_raw")

# Initial filtering-filter missing, minor allele freq, Writing Hardy-Weinberg report (founders only) to firstfilter.hwe 
# --double-id causes both family and within-family IDs to be set to the sample ID.
# --allow-extra-chr allows Nonstandard chromosome IDs
# --geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
# --not-chr is the reverse of --chr: variants on listed chromosome(s) are excluded. 
# --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01)
# # justification here: https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.3
# "If the data will be used for a genome-wide association study (GWAS), usually a rather high MAF threshold is applied, as it requires very strong statistical power to make meaningful statements about very rare alleles. "
# --snps-only excludes all variants with one or more multi-character allele codes
# --hardy writes a list of genotype counts and Hardy-Weinberg equilibrium exact test statistics to plink.hwe

system("./plink --file pearlmusselfiles --double-id --allow-extra-chr --geno 0.2 --mind 0.4 --not-chr NC_043836.1 --maf 0.02 --snps-only --hardy --out firstfilter")
# Total genotyping rate is 0.959533.
## make another ped file with recode to get a .raw file
# system("./plink --file firstfilter --allow-extra-chr --recode A --out firstfilter_raw")


# Now LD prune to get unlinked snps for analysis of population structure
# --indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
# Produces a pruned subset of markers that are in approximate linkage equilibrium with each other, writing the IDs to plink.prune.in (and the IDs of all excluded variants to plink.prune.out). They are currently based on correlations between genotype allele counts; phase is not considered. (Results may be slightly different from PLINK 1.07, due to a minor bugfix in the r2 computation when missing data is present, and more systematic handling of multicollinearity.) Output files are valid input for --extract/--exclude in a future PLINK run.
## below is window size in SNPs (50), the number of SNPs to shift the window at each step (5), the VIF threshold. 0.5 is the R2 threshold for the pairwise SNP-SNP metric (i.e. remove one of a pair of SNPs if the LD is greater than 0.5)
## a. make the prune file:
system("./plink --file pearlmusselfiles --allow-extra-chr --indep-pairwise 50 5 0.5 --out linkage")
#2030 of 5486 variants identified to be removed.
5486-2030 # 3456 remaining
## b. prune the SNPS
system("./plink --file pearlmusselfiles --exclude linkage.prune.out --out LDpruned --make-bed --allow-extra-chr")  
# Total genotyping rate is 0.960745
# 57611 variants and 297 samples pass filters and QC
## now make another .ped file so that you can look at your data
system("./plink --bfile LDpruned --recode tab --out LDpruned_ped --allow-extra-chr")
## make another ped file with recode A so that we get a .raw file
system("./plink --file LDpruned_ped --allow-extra-chr --recode A --out LDpruned_recodeA")

### Heterozygosity
system("./plink --file LDpruned_ped  --out ld_pruned_het --allow-extra-chr --het --recode vcf") #--missing 


## Removing certain populations was tested but that script is removed because it reveals those population's names
