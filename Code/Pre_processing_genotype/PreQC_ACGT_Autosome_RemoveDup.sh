#!/bin/bash

# Pre-processing for genotype data
# PreQC-ACGT,Autosome,RemoveDuplicate for 17 datasets

#1 Keep acgt only for ALL datasets
plink --bfile PITT_3157ppl_925K_hg19 --snps-only just-acgt --make-bed --out PITT_3157ppl_925K_hg19_ACGT
#925265 out of 925517 variants loaded from .bim file.
#3157 people (1310 males, 1847 females) loaded from .fam.
plink --bfile UC_FEMMES_249ppl_593K_hg19_ATGC --snps-only just-acgt --make-bed --out UC_FEMMES_249ppl_593K_hg19_acgt
#593392 out of 593440 variants loaded from .bim file.
#248 people (0 males, 0 females, 248 ambiguous) loaded from .fam.
plink --bfile GHPAFF_3ppl_992K_hg19_ATGC --snps-only just-acgt --make-bed --out GHPAFF_3ppl_992K_hg19_acgt
#991378 out of 992559 variants loaded from .bim file.
#3 people (1 male, 2 females) loaded from .fam.
plink --bfile ADAPT_2784ppl_1M_hg19_ATGC --snps-only just-acgt --make-bed --out ADAPT_2784ppl_1M_hg19_acgt
#1062021 variants loaded from .bim file.
#2784 people (1030 males, 1754 females) loaded from .fam.
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_No3Allele_MergeEqual --snps-only just-acgt --make-bed --out Brazilian_837ppl_1M_hg19_acgt
#1715956 variants loaded from .bim file.
#837 people (280 males, 549 females, 8 ambiguous) loaded from .fam.
plink --bfile CHP_1022ppl_114K_hg19_ATGC --snps-only just-acgt --make-bed --out CHP_1022ppl_114K_hg19_acgt
#114359 out of 114495 variants loaded from .bim file.
#1022 people (333 males, 592 females, 97 ambiguous) loaded from .fam.
plink --bfile CV_hg19_UpdateID_UpdateSex --snps-only just-acgt --make-bed --out CV_697ppl_976K_hg19_acgt
#976939 variants loaded from .bim file.
#697 people (285 males, 412 females) loaded from .fam.
plink --bfile Euro180_176ppl_314K_hg19_ATGC --snps-only just-acgt --make-bed --out Euro180_176ppl_314K_hg19_acgt
#314383 variants loaded from .bim file.
#176 people (0 males, 176 females) loaded from .fam.
plink --bfile FEMMES_256ppl_518K_hg19_ATGC --snps-only just-acgt --make-bed --out FEMMES_256ppl_518K_hg19_acgt
#518476 variants loaded from .bim file.
#255 people (2 males, 253 females) loaded from .fam.
plink --bfile IUPUI_notchr0_ACGT_UpdateName_No3Allele_MergeEqual --snps-only just-acgt --make-bed --out IUPUI_1417ppl_1M_hg19_acgt
#1715941 variants loaded from .bim file.
#1417 people (439 males, 978 females) loaded from .fam.
plink --bfile Lebanese_notchr0_ACGT_UpdateName_No3Allele_MergeEqual --snps-only just-acgt --make-bed --out Lebanese_217ppl_1M_hg19_acgt
#1715941 variants loaded from .bim file.
#217 people (96 males, 121 females) loaded from .fam.
plink --bfile PITT_3157ppl_925K_hg19_ACGT --snps-only just-acgt --make-bed --out PITT_3157ppl_925K_hg19_acgt
#925265 variants loaded from .bim file.
#3157 people (1310 males, 1847 females) loaded from .fam.
plink --bfile SA_231ppl_599K_hg19_ATGC --snps-only just-acgt --make-bed --out SA_231ppl_599K_hg19_acgt
#599830 variants loaded from .bim file.
#231 people (85 males, 146 females) loaded from .fam.
plink --bfile TD2015_199ppl_1M_hg19_ATGC --snps-only just-acgt --make-bed --out TD2015_199ppl_1M_hg19_acgt
#1034392 variants loaded from .bim file.
#199 people (32 males, 167 females) loaded from .fam.
plink --bfile TD2016_notchr0_ACGT_UpdateName_No3Allele_MergeEqual --snps-only just-acgt --make-bed --out TD2016_190ppl_1M_hg19_acgt
#1715941 variants loaded from .bim file.
#190 people (44 males, 146 females) loaded from .fam.
plink --bfile UIUC2013_116ppl_959K_hg19_ATGC --snps-only just-acgt --make-bed --out UIUC2013_116ppl_959K_hg19_acgt
#959801 variants loaded from .bim file.
#116 people (34 males, 82 females) loaded from .fam.
plink --bfile UIUC2014_168ppl_705K_hg19_ATGC --snps-only just-acgt --make-bed --out UIUC2014_168ppl_705K_hg19_acgt
#705914 variants loaded from .bim file.
#168 people (76 males, 92 females) loaded from .fam.
plink --bfile FSUHS_geno0.01_mind0.01_acgt --snps-only just-acgt --make-bed --out FSUHS_57ppl_637K_hg19_acgt
#637397 variants loaded from .bim file.
#57 people (57 males, 0 females) loaded from .fam.

#2 Get autosome only for every dataset
# get all file name
ls | awk -F"." '{print $1}' | uniq > ALL
# plink to extract autosome only for all datasets
for i in $(cat ALL); do plink --bfile $i --autosome --make-bed --out ../3_autosome/${i%.txt}_autosome; done

#3 Find rsids with multiple positions and update
# get all filename
ls | awk -F"." '{print $1}' | uniq > ALL
# get all filename except 4 MEGA chips data in ALL, name it ALL-4
# get rsid and their position for file in ALL-4
for i in $(cat ALL-4); do awk 'BEGIN{OFS="\t"} {print $2,$1,$4}' ${i%.txt}.bim > ${i%.txt}.bim.rsid; done
# get unique position for all rsid
cat *.rsid | awk '!a[$0]++' > Union13_rsid
wc -l Union13_rsid
#2075815 Union13_rsid
# find out if there's any rsid has multiple positions, unique rsids
awk 'BEGIN{OFS="\t"} {print $1}' Union13_rsid | awk '!a[$0]++' | wc -l
#2075706
# about 109 duplicate
awk 'FNR==NR{a[$1]++;next} a[$1]>1' Union13_rsid Union13_rsid | sort -k1,1 > Union13_rsid_dup
# get all duplicate rsid and find hg19 position in ucsc
awk 'BEGIN{OFS="\t"} {print $1}' Union13_rsid_dup | uniq > Union13_rsid_duprsid
# submit to ucsc to get hg19 position of these 109 rsid
#6 of the 109 given identifiers have no match in table snp151, field name. 
#6 missing identifier(s): rs1126793 rs2001943 rs2862633 rs35065085 rs35862547 rs9505819
# for these 6 position can't find rsid in ucsc, get the position in the bim file, 
for i in $(cat union13_missedindb151_6); do grep -w $i *.bim >> miss6_inall; done
#keep the position with more position in datasets
#rs1126793	234526681
#rs2001943	37205505
#rs2862633	176527070
#rs9505819	9200276	
# 2 out of 6 can't decide
#CV_697ppl_976K_hg19_acgt_autosome.bim:7	rs35065085	78.49	57732524	G	A
#PITT_3157ppl_925K_hg19_acgt_autosome.bim:22	rs35065085	0	17246286	G	A
#CV_697ppl_976K_hg19_acgt_autosome.bim:7	rs35862547	87.4411	74837045	A	G
#PITT_3157ppl_925K_hg19_acgt_autosome.bim:7	rs35862547	0	74409382	A	G
# cause in CV, for snp can't find position in ucsc, I used liftover, not sure about if PITT is right. So just keep them.

### re-format to get update file
echo "ucsc<-read.table("union13", header = F, sep="\t")
colnames(ucsc)[1:4]<-c("chr","start","end","rsid")
dup<-read.table("Union13_rsid_dup", header=F, sep="\t")
colnames(dup)[1:3]<-c("rsid","CHR","POS")
MERGE<-merge(ucsc,dup, by=c("rsid"))
MERGE$match<-0
MERGE$match[MERGE$end == MERGE$POS]<-1
#how many match?
table(MERGE$match) #102 MATCH
#re organize the new position
new_ucsc<-ucsc[c("rsid","end")]
write.table(new_ucsc,"Union_UpdatePosition_103", sep="\t",row.names=F, col.names=F, quote=F)" >> UpdatePosition.R
Rscript UpdatePosition.R

#Combine these ucsc updated with the 4 above make Update_Union13 and update with these position for all datasets
for i in $(cat ALL-4); do plink --bfile $i --update-map Update_Union13 --make-bed --out ../4_updatemap/${i%.txt}_updatemap; done
#ADAPT_2784ppl_1M_hg19_acgt_autosome
#--update-map: 78 values updated, 29 variant IDs not present.
#CHP_1022ppl_114K_hg19_acgt_autosome
#--update-map: 2 values updated, 105 variant IDs not present.
#CV_697ppl_976K_hg19_acgt_autosome
#--update-map: 59 values updated, 48 variant IDs not present.
#Euro180_176ppl_314K_hg19_acgt_autosome
#--update-map: 7 values updated, 100 variant IDs not present.
#FEMMES_256ppl_518K_hg19_acgt_autosome
#--update-map: 28 values updated, 79 variant IDs not present.
#FSUHS_57ppl_637K_hg19_acgt_autosome
#--update-map: 6 values updated, 101 variant IDs not present.
#GHPAFF_3ppl_992K_hg19_acgt_autosome
#--update-map: 98 values updated, 9 variant IDs not present.
#PITT_3157ppl_925K_hg19_acgt_autosome
#--update-map: 51 values updated, 56 variant IDs not present.
#SA_231ppl_599K_hg19_acgt_autosome
#--update-map: 36 values updated, 71 variant IDs not present.
#TD2015_199ppl_1M_hg19_acgt_autosome
#--update-map: 70 values updated, 37 variant IDs not present.
#UC_FEMMES_249ppl_593K_hg19_acgt_autosome
#--update-map: 34 values updated, 73 variant IDs not present.
#UIUC2013_116ppl_959K_hg19_acgt_autosome
#--update-map: 70 values updated, 37 variant IDs not present.
#UIUC2014_168ppl_705K_hg19_acgt_autosome
#--update-map: 47 values updated, 60 variant IDs not present.

#4 Change rsid into chr:pos format with --update-name.
# get all filename
ls | awk -F"." '{print $1}' | uniq > ALL
# make update file
for i in $(cat ALL); do awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' ${i%.txt}.bim > ${i%.txt}_updatename; done
# update name
for i in $(cat ALL); do plink --bfile $i --update-name ${i%.txt}_updatename --make-bed --out ../5_updatename/${i%.txt}_updatename; done

#5 Merge those variants with same position to deal with duplicate.
# get all filename
ls | grep 'updatename' | awk -F"." '{print $1}' | uniq > ALL
# find out who is duplicated
for i in $(cat ALL); do awk -F'\t' 'FNR==NR{x[$1,$4]++; next} x[$1,$4]>1' ${i%.txt}.bim ${i%.txt}.bim | sort -k1n,1 -k4n,4 > ${i%.txt}_dup; done
# deal with the datasets with dups
# keep name of datasets  with dups in ALL_dup, 9/13 need merge dups
cat ALL_dup
#ADAPT_2784ppl_1M_hg19_acgt_autosome_updatemap_updatename
#CV_697ppl_976K_hg19_acgt_autosome_updatemap_updatename
#FSUHS_57ppl_637K_hg19_acgt_autosome_updatemap_updatename
#GHPAFF_3ppl_992K_hg19_acgt_autosome_updatemap_updatename
#SA_231ppl_599K_hg19_acgt_autosome_updatemap_updatename
#TD2015_199ppl_1M_hg19_acgt_autosome_updatemap_updatename
#UC_FEMMES_249ppl_593K_hg19_acgt_autosome_updatemap_updatename
#UIUC2013_116ppl_959K_hg19_acgt_autosome_updatemap_updatename
#UIUC2014_168ppl_705K_hg19_acgt_autosome_updatemap_updatename

# use --merge-equal-pos to deal with the duplicate, first merge
for i in $(cat ALL_dup); do plink --bfile $i --bmerge $i --merge-equal-pos --make-bed --out ${i%.txt}_mergeEqualPos; done
#ADAPT_2784ppl_1M_hg19_acgt_autosome_updatemap_updatename
#Error: 1 variant with 3+ alleles present.
#FSUHS_57ppl_637K_hg19_acgt_autosome_updatemap_updatename
#Error: 17 variants with 3+ alleles present.
#SA_231ppl_599K_hg19_acgt_autosome_updatemap_updatename
#Error: 1 variant with 3+ alleles present.
#TD2015_199ppl_1M_hg19_acgt_autosome_updatemap_updatename
#Error: 3 variants with 3+ alleles present.
#UC_FEMMES_249ppl_593K_hg19_acgt_autosome_updatemap_updatename
#Error: 1 variant with 3+ alleles present.
#UIUC2014_168ppl_705K_hg19_acgt_autosome_updatemap_updatename
#Error: 1 variant with 3+ alleles present.
# 3/9 successfully merged

# keep name of datasets with error in ALL_dup2, 6/9 need remove triallele and re-merge
cat ALL_dup2
#ADAPT_2784ppl_1M_hg19_acgt_autosome_updatemap_updatename
#FSUHS_57ppl_637K_hg19_acgt_autosome_updatemap_updatename
#SA_231ppl_599K_hg19_acgt_autosome_updatemap_updatename
#TD2015_199ppl_1M_hg19_acgt_autosome_updatemap_updatename
#UC_FEMMES_249ppl_593K_hg19_acgt_autosome_updatemap_updatename
#UIUC2014_168ppl_705K_hg19_acgt_autosome_updatemap_updatename

# exclude triallele
for i in $(cat ALL_dup2); do plink --bfile $i --exclude ${i%.txt}_mergeEqualPos-merge.missnp --make-bed --out ${i%.txt}_No3Allele; done
# make the name list of these 6 datasets
ls | grep 'No3Allele' | awk -F"." '{print $1}' | uniq > ALL_dup3
# second merge
for i in $(cat ALL_dup3); do plink --bfile $i --bmerge $i --merge-equal-pos --make-bed --out ../6_mergeEqualPos/${i%.txt}_mergeEqualPos; done
