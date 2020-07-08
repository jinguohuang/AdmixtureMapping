#!/bin/bash

# Ancestry estimation
## QC after merging and LD based SNP pruning


for name in $(cat ALL) 
do
	# QC for merged data 
	echo "QC for ${name}"
	plink --bfile ${name} --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out ${name}_geno_mind_maf
	# LD pruning after QC
	echo "LD pruning for ${name}"
	plink --bfile ${name}_geno_mind_maf --indep-pairwise 100 10 0.1 --out ${name}_geno_mind_maf
	# Extract pruned SNP
	echo "Extracting pruned SNP for ${name}"
	plink --bfile ${name}_geno_mind_maf --extract ${name}_geno_mind_maf.prune.in --make-bed --out ../5_LDpruning/${name}_geno_mind_maf_pruned
done

