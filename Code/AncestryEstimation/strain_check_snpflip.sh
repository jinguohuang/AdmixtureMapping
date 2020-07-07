#!/bin/bash

# Ancestry estimation
## Strain check with snpflip

# Reference data preparation: 1KG dataset in fasta format
# Download data from 1KG website
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
# Decompress the file
gunzip human_g1k_v37.fasta.gz


# snpflip software installation (https://github.com/biocore-ntnu/snpflip)
pip install snpflip

# Our data have been through QC (--geno 0.1 --mind 0.1 remove palindromic and --maf 0.05)
# List filename of our dataset in a file so we can use them later with command
ls | rev | cut -f 2- -d '.' | rev | uniq > ALL
# Check what's inside the file
cat ALL
'''
ADAPT_preQC_geno_mind_NoPalindromic_maf0.05
Brazilian_preQC_geno_mind_NoPalindromic_maf0.05
CHP_preQC_geno_mind_NoPalindromic_maf0.05
CV_preQC_geno_mind_NoPalindromic_maf0.05
Euro180_preQC_geno_mind_NoPalindromic_maf0.05
FEMMES_preQC_geno_mind_NoPalindromic_maf0.05
FSUHS_preQC_geno_mind_NoPalindromic_maf0.05
GENSET_v4_NoPalindromic_maf0.05
GENSET_v5_NoPalindromic_maf0.05
GHPAFF_preQC_geno_mind_NoPalindromic_maf0.05
IUPUI_preQC_geno_mind_NoPalindromic_maf0.05
Lebanese_preQC_geno_mind_NoPalindromic_maf0.05
PITT_preQC_geno_mind_NoPalindromic_maf0.05
SA_preQC_geno_mind_NoPalindromic_maf0.05
TD2015_preQC_geno_mind_NoPalindromic_maf0.05
TD2016_preQC_geno_mind_NoPalindromic_maf0.05
UIUC2013_preQC_geno_mind_NoPalindromic_maf0.05
UIUC2014_preQC_geno_mind_NoPalindromic_maf0.05
'''

# Run snpflip to check strain of our datasets against the reference data
for i in $(cat ALL)
do 
	snpflip --fasta-genome ../0_ReferenceData/1000G_hg19_fasta/human_g1k_v37.fasta \
	--bim-file ${i%.txt}.bim \
	--output-prefix $i
done

# Run plink to flip reverse SNP and remove ambiguous SNP.
for i in $(cat ALL) 
do 
	plink --bfile $i \
	--flip ${i%.txt}.reverse \
	--exclude ${i%.txt}.ambiguous \
	--make-bed --out ${i%.txt}_Snpflip
done



