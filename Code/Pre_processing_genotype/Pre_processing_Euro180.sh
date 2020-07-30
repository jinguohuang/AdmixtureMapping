#!/bin/bash

# Pre-processing for genotype data - Euro180
## Euro180 was hg17 and need updated with hg19 positions

#1 Get hg19 position for all variants
# Get rsid and submit to get UCSC hg19 position
awk 'BEGIN{OFS="\t"}{print $2}' Euro180_176ppl_317K_hg19_ATGC.bim > Euro180_rsid
#3413 of the 317503 given identifiers have no match in table snp151Common
#3114 of the 3413 given identifiers have no match in table snp151

#2 Prepare position update list and remove list
Rscript rsid_update_Euro180.R

#3 Update Plink file
plink --bfile Euro180_176ppl_317K_hg19_ATGC --update-map Euro180_update_map --exclude Euro180_remove_rsid --make-bed --out Euro180_176ppl_314K_hg19_ATGC
#317503 variants loaded from .bim file.
#176 people (0 males, 176 females) loaded from .fam.
#--update-map: 314383 values updated.
#Warning: Base-pair positions are now unsorted!
#--exclude: 314383 variants remaining.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.991986.
#314383 variants and 176 people pass filters and QC.
