#!/bin/bash

# Pre-processing for genotype data
## Process IUPUI, TD2016, Brazilian, Lebanese (MEGA chips)

#1 Remove chromosome 0
## TD2016
plink --bfile NewTD --not-chr 0 --make-bed --out TD2016_notchr0
#1768351 out of 1779819 variants loaded from .bim file.
#190 people (44 males, 146 females) loaded from .fam.
#Warning: 20328 het. haploid genotypes present (see TD2016_notchr0.hh ); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.990241.
#1768351 variants and 190 people pass filters and QC.

## IUPUI
plink --bfile PennStateGenotypes --not-chr 0 --make-bed --out IUPUI_notchr0
#1768351 out of 1779819 variants loaded from .bim file.
#1417 people (439 males, 978 females) loaded from .fam.
#Warning: 173018 het. haploid genotypes present (see IUPUI_notchr0.hh ); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.995429.
#1768351 variants and 1417 people pass filters and QC.

## Lebanese
plink --bfile PSU_Lebanese --not-chr 0 --make-bed --out Lebanese_notchr0
#1768351 out of 1779819 variants loaded from .bim file.
#217 people (96 males, 121 females) loaded from .fam.
#Warning: 34961 het. haploid genotypes present (see Lebanese_notchr0.hh ); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.998616.
#1768351 variants and 217 people pass filters and QC.

## Brazilian
plink --bfile PSU_Brazilians --not-chr 0 --make-bed --out Brazilian_notchr0
#1768351 out of 1779819 variants loaded from .bim file.
#837 people (0 males, 0 females, 837 ambiguous) loaded from .fam.
#Ambiguous sex IDs written to Brazilian_notchr0.nosex .
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.985697.
#1768351 variants and 837 people pass filters and QC.

#2 Brazilian update sex, update id
#Update Brazilian Sex
plink --bfile PSU_Brazilians_SexUpdated --update-sex Braz_SexUpdate.txt --make-bed --out Brazilians_SexUpdated          
#update id from 94594 to 64594
plink --bfile Brazilians_SexUpdated --update-ids Braz_IDUpdate.txt --make-bed --out Brazilians_SexUpdated_IDUpdated
# --update-ids: 1 person updated.

#3 Remove deletion or insertion
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0 --snps-only just-acgt --make-bed --out Brazilians_SexUpdated_IDUpdated_notchr0_ACGT
#1742804 out of 1768351 variants loaded from .bim file.
plink --bfile Lebanese_notchr0 --snps-only just-acgt --make-bed --out Lebanese_notchr0_ACGT
#1742804 out of 1768351 variants loaded from .bim file.
plink --bfile IUPUI_notchr0 --snps-only just-acgt --make-bed --out IUPUI_notchr0_ACGT
#1742804 out of 1768351 variants loaded from .bim file.
plink --bfile TD2016_notchr0 --snps-only just-acgt --make-bed --out TD2016_notchr0_ACGT
#1742804 out of 1768351 variants loaded from .bim file.

#4 Change rsID into chr:pos format and remove duplicated position
#find out who is dup
awk -F'\t' 'FNR==NR { x[$1,$4]++; next } x[$1,$4] > 1' Lebanese_notchr0_ACGT.bim Lebanese_notchr0_ACGT.bim | sort -k1n,1 -k4n,4 > Lebanese_dup
wc -l Lebanese_dup
#52603 Lebanese_dup
# use --merge-equal-pos to deal with the duplicate

#Lebanese
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' Lebanese_notchr0_ACGT.bim > Lebanese_UpdateName.txt
plink --bfile Lebanese_notchr0_ACGT --update-name Lebanese_UpdateName.txt --make-bed --out Lebanese_notchr0_ACGT_UpdateName
#--update-name: 1742804 values updated.
plink --bfile Lebanese_notchr0_ACGT_UpdateName --bmerge Lebanese_notchr0_ACGT_UpdateName --merge-equal-pos --make-bed --out Lebanese_notchr0_ACGT_UpdateName_mergeEqualPos
#Error: 515 variants with 3+ alleles present.
# remove trialleles
plink --bfile Lebanese_notchr0_ACGT_UpdateName --exclude RemoveTheseTriallele.txt --make-bed --out Lebanese_notchr0_ACGT_UpdateName_No3Allele
#1742804 variants loaded from .bim file.
#--exclude: 1741769 variants remaining.
# merge again
plink --bfile Lebanese_notchr0_ACGT_UpdateName_No3Allele --bmerge Lebanese_notchr0_ACGT_UpdateName_No3Allele --merge-equal-pos --make-bed --out Lebanese_notchr0_ACGT_UpdateName_No3Allele_MergeEqual

# other 3 datasets
# update name list
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' Brazilians_SexUpdated_IDUpdated_notchr0_ACGT.bim > Braz_UpdateName.txt
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' IUPUI_notchr0_ACGT.bim > IUPUI_UpdateName.txt
awk 'BEGIN{OFS="\t"}{print $2,$1":"$4}' TD2016_notchr0_ACGT.bim > TD2016_UpdateName.txt
# update name
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0_ACGT --update-name Braz_UpdateName.txt --make-bed --out Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName
#--update-name: 1742804 values updated.
plink --bfile IUPUI_notchr0_ACGT --update-name IUPUI_UpdateName.txt --make-bed --out IUPUI_notchr0_ACGT_UpdateName
#--update-name: 1742804 values updated.
plink --bfile TD2016_notchr0_ACGT --update-name TD2016_UpdateName.txt --make-bed --out TD2016_notchr0_ACGT_UpdateName
#--update-name: 1742804 values updated.
# first merge to get tri allele list
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName --bmerge Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName --merge-equal-pos --make-bed --out Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_MergeEqual
#Error: 500 variants with 3+ alleles present.
mv Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_MergeEqual-merge.missnp RemoveTheseTriallele_Braz.txt
plink --bfile IUPUI_notchr0_ACGT_UpdateName --bmerge IUPUI_notchr0_ACGT_UpdateName --merge-equal-pos --make-bed --out IUPUI_notchr0_ACGT_UpdateName_MergeEqual
#Error: 515 variants with 3+ alleles present.
mv IUPUI_notchr0_ACGT_UpdateName_MergeEqual-merge.missnp RemoveTheseTriallele_IUPUI.txt
plink --bfile TD2016_notchr0_ACGT_UpdateName --bmerge TD2016_notchr0_ACGT_UpdateName --merge-equal-pos --make-bed --out TD2016_notchr0_ACGT_UpdateName_MergeEqual
#Error: 515 variants with 3+ alleles present.
mv TD2016_notchr0_ACGT_UpdateName_MergeEqual-merge.missnp RemoveTheseTriallele_TD.txt

#remove those tri allele
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName --exclude RemoveTheseTriallele_Braz.txt --make-bed --out Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_No3Allele
#1742804 variants loaded from .bim file.
#--exclude: 1741799 variants remaining.
plink --bfile IUPUI_notchr0_ACGT_UpdateName --exclude RemoveTheseTriallele_IUPUI.txt --make-bed --out IUPUI_notchr0_ACGT_UpdateName_No3Allele
#1742804 variants loaded from .bim file.
#--exclude: 1741769 variants remaining.
plink --bfile TD2016_notchr0_ACGT_UpdateName --exclude RemoveTheseTriallele_TD.txt --make-bed --out TD2016_notchr0_ACGT_UpdateName_No3Allele
#1742804 variants loaded from .bim file.
#--exclude: 1741769 variants remaining.

# second merge to remove the duplicate
plink --bfile Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_No3Allele --bmerge Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_No3Allele --merge-equal-pos --make-bed --out Brazilians_SexUpdated_IDUpdated_notchr0_ACGT_UpdateName_No3Allele_MergeEqual
plink --bfile IUPUI_notchr0_ACGT_UpdateName_No3Allele --bmerge IUPUI_notchr0_ACGT_UpdateName_No3Allele --merge-equal-pos --make-bed --out IUPUI_notchr0_ACGT_UpdateName_No3Allele_MergeEqual
plink --bfile TD2016_notchr0_ACGT_UpdateName_No3Allele --bmerge TD2016_notchr0_ACGT_UpdateName_No3Allele --merge-equal-pos --make-bed --out TD2016_notchr0_ACGT_UpdateName_No3Allele_MergeEqual
