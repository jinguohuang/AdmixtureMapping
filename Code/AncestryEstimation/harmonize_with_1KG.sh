#!/bin/bash

# Ancestry estimation
## Harmonize datasets with reference data using Genotype Harmonizer

# i. Filter out ambiguous SNPs
# ii. Update ids to match 1000G
# iii. Update reference allele
# iv. Merge harmonized chromosome files

# Reference data preparation: 1KG dataset in vcf format
# Run script for downloading 1KG vcf
chmod +x 1KG_vcf_downloader.sh
./1KG_vcf_downloader.sh

# Genotype Harmonizer software installation
wget https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.zip
unzip GenotypeHarmonizer-1.4.23-dist.zip

# Our data have been though snpflip in the last step
# Prepare a file that list all filenames without suffix
'''
ADAPT_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
Brazilian_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
CHP_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
CV_preQC_geno_mind_NoPalindromic_maf0.05_exclude_snpflip
Euro180_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
FEMMES_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
FSUHS_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
GENSET_v4_NoPalindromic_maf0.05_snpflip
GENSET_v5_NoPalindromic_maf0.05_snpflip
GHPAFF_preQC_geno_mind_NoPalindromic_maf0.05_exclude_snpflip
IUPUI_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
Lebanese_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
PITT_preQC_geno_mind_NoPalindromic_maf0.05_exclude_snpflip
SA_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
TD2015_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
TD2016_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
UIUC2013_preQC_geno_mind_NoPalindromic_maf0.05_exclude_snpflip
UIUC2014_preQC_geno_mind_NoPalindromic_maf0.05_snpflip
'''


#for each file, split them by chromosome and run the harmonizer
thisdir=$(pwd)
for name in $(cat ALL) 
do
	#name= `echo "$file"` 
	mkdir $name
	mkdir $name"_temp"
	echo "Spliting file $name by chromosome"
	
	for chr in {1..22}
	do
		plink --bfile $name --chr $chr --make-bed --out ${name}_temp/${name}_chr${chr}
	done

	#Create one job for each file for each chromosome
	#This step will take a while
	for i in {1..22}
	do
		echo "Writing job for $name"
		echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l pmem=8gb
#PBS -A tll30_a_g_bc_default
#PBS -j oe
cd ${thisdir}
java -Xmx2048m -jar  GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar --input ${name}_temp/${name}_chr${i} \
--inputType PLINK_BED \
--ref ~/scratch/1KG_vcf/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
--refType VCF \
--ambiguousSnpFilter \
--update-id \
--min-ld 0.3 \
--mafAlign 0.1 \
--min-variants 3 \
--variants 100 \
--update-reference-allele \
--keep \
--debug \
--output ${name}/${name}_chr${i}_harmonized" >> job_${name}_chr${i}.pbs

		echo "Submitting job_${name}_chr${i}.pbs"
		qsub job_${name}_chr${i}.pbs
		echo "Waiting 5s for next chromosome..."
		sleep 5s 

	done
	echo "Waiting 15m for next file..."
	sleep 15m
done

# Merging all chromosome together for each dataset after harmonizing.
for name in $(cat ALL) 
do
	echo "Merging file $name by chromosome"
	for chr in {2..22}
	do
		echo "${name}/${name}_chr${chr}_harmonized" >> ${name}_chr_MergeList.txt
	done
	plink --bfile ${name}/${name}_chr1_harmonized \
	--merge-list ${name}_chr_MergeList.txt \
	--make-bed --out ${name}/${name}_HarmonizedTo1000G
done


