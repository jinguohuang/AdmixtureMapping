#!/bin/bash

# Pre-processing for genotype data - GENSET

#1 Convert raw data from vcf to plink format
plink --vcf GENSET_5_17_2018.vcf --const-fid GENSET --biallelic-only strict list --make-bed --out GENESET
#1487373 variants loaded from .bim file.
#248 people (0 males, 0 females, 248 ambiguous) loaded from .fam.
#Total genotyping rate is 0.40556.
#1487373 variants and 248 people pass filters and QC.

#2 Identify SNPs from v4 and v5 chips of 23andme
##2.1 Check the SNP missingness frequency
plink --bfile GENSET --missing --out GENSET
#The .lmiss files from plink have a variable number of whitespaces as delimiter
#convert variable whitespace to a single whitespace.
sed ' s/ \+ / /g' GENSET.lmiss > GENSET_space.lmiss

##2.2 Plot histogram for GENSET missingness and split the dataset to make SNP list for v4 v5
echo"
library(readr)
library(ggplot2)
lmiss<-read.table("GENSET_space.lmiss",sep = " ", header=T)
#Make a plot
p<-ggplot(lmiss, aes(F_MISS)) +
  geom_histogram(binwidth=0.01)+
  theme_bw()+
  ggtitle("SNP missingness in GENSET")
ggsave(plot = p, filename = "GENSET_SNPmissingness.jpg", device = "jpeg", width = 6, height = 4, units="in")
# 
MissingInAll <- lmiss[which(lmiss$F_MISS >= 0.75),]
MissingIn155 <- lmiss[which(lmiss$F_MISS >= 0.5 & lmiss$F_MISS <= 0.75),]
MissingIn93 <- lmiss[which(lmiss$F_MISS >= 0.25 & lmiss$F_MISS <= 0.5),]
MissingInNone <- lmiss[which(lmiss$F_MISS <= 0.25),]
write.table(MissingInAll$SNP, "MissingInAll.txt", row.names = F, col.names = F, quote = F) #382149
write.table(MissingIn155$SNP, "MissingIn155.txt", row.names = F, col.names = F, quote = F) #506985
write.table(MissingIn93$SNP, "MissingIn93.txt", row.names = F, col.names = F, quote = F) #479484
write.table(MissingInNone$SNP, "MissingInNone.txt", row.names = F, col.names = F, quote = F) #112917
" >> GENSET_split.R
## run
Rscript GENSET_split.R

##2.3 make SNP lists for splitting the dataset into the v4 and v5 arrays
cat MissingIn93.txt MissingInNone.txt > GENSET_v4_SNPs.txt
cat MissingIn155.txt MissingInNone.txt > GENSET_v5_SNPs.txt
plink --bfile GENSET --extract GENSET_v4_SNPs.txt --mind 0.1 --geno 0.1 --make-bed --out GENSET_v4
#589739 variants and 155 people pass filters and QC.
plink --bfile GENSET --extract GENSET_v5_SNPs.txt --mind 0.1 --geno 0.1 --make-bed --out GENSET_v5
#617490 variants and 93 people pass filters and QC.

#3 Update individual IDs
plink --bfile GENSET_v4 --update-ids GenSet_UpdateIID.txt --make-bed --out ../3_UpdateIDs/GENSET_v4_UpdateID
#--update-ids: 155 people updated, 96 IDs not present.
plink --bfile GENSET_v5 --update-ids GenSet_UpdateIID.txt --make-bed --out ../3_UpdateIDs/GENSET_v5_UpdateID
#--update-ids: 93 people updated, 158 IDs not present.

#4 Update sex
##4.1 Remove MT chromosome for sex calculation
plink --bfile GENSET_v4_UpdateID --chr 1-24 --make-bed --out ../4_RemoveMT/GENSET_v4_UpdateID_MTRemoved
plink --bfile GENSET_v5_UpdateID --chr 1-24 --make-bed --out ../4_RemoveMT/GENSET_v5_UpdateID_MTRemoved
##4.2 Calculate sex
# split psudoautosome 
plink --bfile GENSET_v4_UpdateID_MTRemoved --split-x hg19 no-fail --make-bed --out ../5_UpdateSex/GENSET_v4_UpdateID_MTRemoved_splitX
plink --bfile GENSET_v5_UpdateID_MTRemoved --split-x hg19 no-fail --make-bed --out ../5_UpdateSex/GENSET_v5_UpdateID_MTRemoved_splitX
# LD pruning
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX --indep-pairwise 50 5 0.5 --out GENSET_v4_UpdateID_MTRemoved_splitX
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX --indep-pairwise 50 5 0.5 --out GENSET_v5_UpdateID_MTRemoved_splitX
# sex check
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX --exclude GENSET_v4_UpdateID_MTRemoved_splitX.prune.out --check-sex 0.3 0.7 --out GENSET_v4_UpdateID_MTRemoved_splitX
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX --exclude GENSET_v5_UpdateID_MTRemoved_splitX.prune.out --check-sex 0.3 0.7 --out GENSET_v5_UpdateID_MTRemoved_splitX
# make update sex list
awk -v OFS="\t" 'NR>1 {print $1,$2,$4}' GENSET_v4_UpdateID_MTRemoved_splitX.sexcheck > GENSET_v4_UpdateSex.txt
awk -v OFS="\t" 'NR>1 {print $1,$2,$4}' GENSET_v5_UpdateID_MTRemoved_splitX.sexcheck > GENSET_v5_UpdateSex.txt
#Update sex with plink.
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX --update-sex GENSET_v4_UpdateSex.txt --make-bed --out GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX --update-sex GENSET_v5_UpdateSex.txt --make-bed --out GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex

#5 Keep acgt only 
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex --snps-only just-acgt --make-bed --out ../6_acgt_Autosome/GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex --snps-only just-acgt --make-bed --out ../6_acgt_Autosome/GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt

#6 Keep autosome only.
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt --autosome --make-bed --out GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt --autosome --make-bed --out GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome

#7 Change rsid into chr:pos format
##7.1 make the update name list first
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome.bim > GENSET_V4_updatename
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome.bim > GENSET_V5_updatename
##7.2 Update name
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome --update-name GENSET_V4_updatename --make-bed --out ../7_UpdateName/GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome --update-name GENSET_V5_updatename --make-bed --out ../7_UpdateName/GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName

#8 Check duplicate position
awk -F"\t" '{print $2}' GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName.bim | sort | uniq | wc -l
#594851
awk -F"\t" '{print $2}' GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName.bim | sort | uniq | wc -l
#564407
# no duplicate

#9 Remove palindromic SNPs
#Get GC/CG/AT/TA position in each dataset 
grep -P $'G\tC'\|$'C\tG'\|$'A\tT'\|$'T\tA' GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName.bim > GENSET_V4_palindromicSNP
grep -P $'G\tC'\|$'C\tG'\|$'A\tT'\|$'T\tA' GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName.bim > GENSET_V5_palindromicSNP
awk -F"\t" '{print $2}' GENSET_V4_palindromicSNP > GENSET_V4_palindromicSNP_rsid
awk -F"\t" '{print $2}' GENSET_V5_palindromicSNP > GENSET_V5_palindromicSNP_rsid
#Remove them
plink --bfile GENSET_v4_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName --exclude GENSET_V4_palindromicSNP_rsid --make-bed --out ../8_RemovePalindromicSNP/GENSET_v4_NoPalindromic
plink --bfile GENSET_v5_UpdateID_MTRemoved_splitX_UpdateSex_acgt_autosome_UpdateName --exclude GENSET_V5_palindromicSNP_rsid --make-bed --out ../8_RemovePalindromicSNP/GENSET_v5_NoPalindromic
# Done