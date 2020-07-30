#!/bin/bash

# Pre-processing for genotype data - CV
## LiftOver from hg18 to hg19 and update individual id
## Update all rsid with hg19 positions, update all non-rsid and rsid that not found match in UCSC with liftOver, remove the rest


# Start with raw file CapeVerde.ped and CapeVerde.map
# Convert raw ped file into plink binary file
plink --file CapeVerde --make-bed --out liftOver/CV_hg18

# get rsid and non-rsid
grep "rs" CV_hg18.bim > CV_hg18.bim.rs
grep -v "rs" CV_hg18.bim > CV_hg18.bim.nonrs

#1 non-rsid position: liftOver from hg18 to hg19
#1.1 convert bim file into bed file so liftOver can process
Rscript bim_to_bed.R CV_hg18.bim.nonrs
#1.2 liftover
# get chain file for liftover first
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
liftOver CV_hg18.bim.nonrs_bimtobed hg18ToHg19.over.chain CV_hg18.bim.nonrs_bimtobed_hg19 CV_hg18.bim.nonrs_bimtobed_unmapped
# check the result files
wc -l CV_hg18.bim.nonrs_bimtobed*
  71153 CV_hg18.bim.nonrs_bimtobed
  71145 CV_hg18.bim.nonrs_bimtobed_hg19
     16 CV_hg18.bim.nonrs_bimtobed_unmapped
# 71145 mapped, 8 unmapped
#1.3 make update list from mapped file
awk 'BEGIN{OFS="\t"}; {print $4,$3}' CV_hg18.bim.nonrs_bimtobed_hg19 > update_nonrs.txt
#1.4 make exclude list from unmapped file
grep -v '^#' CV_hg18.bim.nonrs_bimtobed_unmapped | awk 'BEGIN{OFS="\t"}; {print $4}' > remove_nonrs.txt

#2 rsid position: update using UCSC hg19 position
#2.1 Get hg19 position for all rsid from UCSC
# get rsid
awk 'BEGIN{OFS="\t"}; {print $2}' CV_hg18.bim.rs > CV_hg18.bim.rsid
# upload to ucsc
# 43464 of the 906322 given identifiers have no match in table snp151Common
# split into 5 files and submit to 151all
split -l 8000 CV_common151_missed CV_common151_missed_
#2.2 combine all CV rsid hg19 pos and Update positions that have same chromosome number
Rscript rsid_update_CV.R

#3 rsid not in UCSC: liftOver positions
# bim to bed
Rscript bim_to_bed.R CV_hg18.bim.rs.removed
# liftOver
liftOver CV_hg18.bim.rs.removed_bimtobed hg18ToHg19.over.chain CV_hg18.bim.rs.removed_bimtobed_hg19 CV_hg18.bim.rs.removed_bimtobed_unmapped
wc -l CV_hg18.bim.rs.removed_bimtobed*
  2010 CV_hg18.bim.rs.removed_bimtobed
  1907 CV_hg18.bim.rs.removed_bimtobed_hg19
   206 CV_hg18.bim.rs.removed_bimtobed_unmapped
# liftOver 1907, 103 fail
# make update list and remove list
awk 'BEGIN{OFS="\t"}; {print $4,$3}' CV_hg18.bim.rs.removed_bimtobed_hg19 > update_rs_removed.txt
grep -v '^#' CV_hg18.bim.rs.removed_bimtobed_unmapped | awk 'BEGIN{OFS="\t"}; {print $4}' > remove_rs_removed.txt

#4 Combine update and remove list
#4.1 update list
cat CV_hg18.bim.rs.update update_rs_removed.txt update_nonrs.txt > update_cvhg19.txt
wc -l update_cvhg19.txt
# 977364 update_cvhg19.txt
#4.2 remove list
cat remove_nonrs.txt remove_rs_removed.txt > remove_cv.txt
wc -l remove_cv.txt
#111 remove_cv.txt

#5 Update genotype files and convert individual ids
#5.1 --snps-only just-acgt remove D/I
plink --bfile CV_hg18 --snps-only just-acgt --update-map update_cvhg19.txt --exclude remove_cv.txt --make-bed --out CV_hg19
#977050 out of 977475 variants loaded from .bim file.
#697 people (0 males, 0 females, 697 ambiguous) loaded from .fam.
#--update-map: 976939 values updated, 425 variant IDs not present.
#Warning: Base-pair positions are now unsorted!
#--exclude: 976939 variants remaining.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.997534.
#976939 variants and 697 people pass filters and QC.
#5.2 Convert ids for individuals
plink --bfile CV_hg19 --update-ids CV_ConvertID_new.txt --make-bed --out CV_hg19_UpdateID
#697 people (0 males, 0 females, 697 ambiguous) loaded from .fam.
#--update-ids: 697 people updated.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
#Total genotyping rate is 0.997534.