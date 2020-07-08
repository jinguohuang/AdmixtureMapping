#!/bin/bash

# Ancestry estimation
## Merge our datasets with 1KG and QC

# This scirpt is for extracting SNP in datasets from 1KG plink file and merging our dataset with 1KG

#first make SNP list 

thisdir=$(pwd)

for name in $(cat ALL) 
do
	# prepare harmonized datasets
	echo "Copying $name harmonized file to this folder" 
	cp ../2_snpflip/${name}/${name}_HarmonizedTo1000G.* .
	# prepare SNP in our dataset for extracting
	echo "Getting SNP from $name"
	awk '{print $2}' ${name}_HarmonizedTo1000G.bim | grep rs > ${name}_SNP
	# extract SNP from 1KG for better merging
	echo "Extracting SNP of $name from 1KG"
	plink --bfile all_phase3 --allow-extra-chr --extract ${name}_SNP --make-bed --out all_phase3_extract_${name}
	# merge our data with 1KG
	echo "Merging $name and 1KG"
	plink --bfile ${name}_HarmonizedTo1000G --bmerge all_phase3_extract_${name} --make-bed --out ${name}_1000G_merge
done

