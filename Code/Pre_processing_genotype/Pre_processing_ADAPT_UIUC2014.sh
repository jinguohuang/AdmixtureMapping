#!/bin/bash

# Pre-processing for genotype data - ADAPT, UIUC2014
## LiftOver from hg18 to hg19 and merge 

#1 Remove D/I in raw file (deletion or insertion), this can also be done with plink --snp-only just-acgt
# Load within folder of 23andme format raw genotype data of each individuals in ADAPT or UIUC2014
# delete all comment lines and all rows contains D/I
for i in 1*.txt; do grep -v '^#' "$i" | grep -v "I\|D" > ../UIUC2014_V3/$i.noDI; done
for i in 1*.txt; do grep -v '^#' "$i" | grep -v "I\|D" > ../ADAPT_V3/$i.noDI; done 
# the above command may take long time in ADAPT, so write a pbs file to run on cluster
echo "for i in 1*.txt
do grep -v '^#' "$i" | grep -v "I\|D" > ../ADAPT_V3/$i.noDI
done" >> removeDI.sh
# make it executable
chmod +x removeDI.sh
# pbs file
echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -j oe
cd $PBS_O_WORKDIR
bash removeDI.sh" >> removeDI.pbs
# submit to cluster
qsub -A open removeDI.pbs

#2 Update 132001, 140603 and 140669 from hg18 to hg19 positions
#2.0 prepare for update, get uniq position list for 167 individuals in UIUC2014
cat 14*.noDI | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk '!a[$0]++' > UIUC2014_167_noDI_uniqPOS
# divide files by rsid and iid, treat them separately
# in UIUC2014
grep "i" 132001.txt.noDI > 132001.txt.noDI.i
grep "rs" 132001.txt.noDI > 132001.txt.noDI.rs
# in ADAPT
grep "i" 140603.txt.noDI > 140603.txt.noDI.i
grep "rs" 140603.txt.noDI > 140603.txt.noDI.rs
grep "i" 140669.txt.noDI > 140669.txt.noDI.i
grep "rs" 140669.txt.noDI > 140669.txt.noDI.rs
# get iid position from 167ppl unique position and update
grep "i" UIUC2014_167_noDI_uniqPOS > UIUC2014_167_noDI_uniqPOS_i

##2.1 Update iid and rsid position using other individuals in UIUC2014
## change in the script for 3 individuals
Rscript rsid_update_ADAPT_UIUC2014.R

##2.2 Combine rsid and iid updated positions
cat 132001.txt.noDI.i.update 140603.txt.noDI.rs.update >> 132001.txt.noDI.hg19update
cat 140603.txt.noDI.i.update 140603.txt.noDI.rs.update >> 140603.txt.noDI.hg19update
cat 140669.txt.noDI.i.update 140669.txt.noDI.rs.update >> 140669.txt.noDI.hg19update
##2.3 Sort them by chromosome 1-22,X,Y,MT, so later plink merge would not warning
sort -Vk2,2 -k3n,3 132001.txt.noDI.hg19update > 132001.txt.noDI.hg19update.sortV
sort -Vk2,2 -k3n,3 140603.txt.noDI.hg19update > 140603.txt.noDI.hg19update.sortV
sort -Vk2,2 -k3n,3 140669.txt.noDI.hg19update > 140669.txt.noDI.hg19update.sortV
grep -v "MT" 132001.txt.noDI.hg19update.sortV > 132001.txt.noDI.hg19update.sortV.noMT
grep "MT" 132001.txt.noDI.hg19update.sortV > 132001.txt.noDI.hg19update.sortV.MT
grep -v "MT" 140603.txt.noDI.hg19update.sortV > 140603.txt.noDI.hg19update.sortV.noMT
grep "MT" 140603.txt.noDI.hg19update.sortV > 140603.txt.noDI.hg19update.sortV.MT
grep -v "MT" 140669.txt.noDI.hg19update.sortV > 140669.txt.noDI.hg19update.sortV.noMT
grep "MT" 140669.txt.noDI.hg19update.sortV > 140669.txt.noDI.hg19update.sortV.MT
cat 132001.txt.noDI.hg19update.sortV.noMT 132001.txt.noDI.hg19update.sortV.MT > 132001.txt.noDI
cat 140603.txt.noDI.hg19update.sortV.noMT 140603.txt.noDI.hg19update.sortV.MT > 140603.txt.noDI
cat 140669.txt.noDI.hg19update.sortV.noMT 140669.txt.noDI.hg19update.sortV.MT > 140669.txt.noDI

#3 Convert from 23andme to Plink format and merge
# ADAPT
##3.1 Create update list for mismatch positions
# Check if there's any duplicate or mismatch of rsid and pos
cat *.noDI | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk '!a[$0]++' > ADAPT_ALL.noDI.uniqPOS
wc -l ADAPT_ALL.noDI.uniqPOS
#1079136 ADAPT_ALL.noDI.uniqPOS
# Check if there's any duplicate/mismatch of rsid and pos
awk 'BEGIN{OFS="\t"}{print $1}' ADAPT_ALL.noDI.uniqPOS | awk '!a[$0]++' | wc -l
#1079126
#20 mismatch, who are they?
awk -F'\t' 'FNR==NR {x[$1]++; next} x[$1]>1' ADAPT_ALL.noDI.uniqPOS ADAPT_ALL.noDI.uniqPOS
#check the position for them in the dbSNP
awk -F'\t' 'FNR==NR {x[$1]++; next} x[$1]>1' ADAPT_ALL.noDI.uniqPOS ADAPT_ALL.noDI.uniqPOS | awk -F'\t' '{print $1}' | sort | uniq -dc
# create update map.
echo "rs10929316	234070897
rs2155163	126524251
rs8056335	28631334
rs14968	19956434
rs5758598	42538462
rs2032678	15508700
rs3906	21724699
rs2032634	21888873
rs2032675	21894172
rs2032615	21930259
" >> UpdateTheseAllele.txt
##3.2 Create remove list for duplicated iid positions
# check how many have duplicated positions
awk -F'\t' 'FNR==NR {x[$2,$3]++; next} x[$2,$3]>1' ADAPT_ALL.noDI.uniqPOS ADAPT_ALL.noDI.uniqPOS > ADAPT_ALL.noDI.dupPOS
wc -l ADAPT_ALL.noDI.dupPOS
#25111 ADAPT_ALL.noDI.dupPOS
grep "i" ADAPT_ALL.noDI.dupPOS | awk -F'\t' '{print $1}' > ADAPT_ALL.noDI.dupPOS.iid
# make the iid list part of remove list
mv ADAPT_ALL.noDI.dupPOS.iid RemoveTheseAllele_V2.txt

##3.3 First merge trial and update the remove list
echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -j oe
cd $PBS_O_WORKDIR
./ADAPT_23andMeToPlink_V2.bat"
chmod u+x ADAPT_MERGE_V2.pbs
qsub -A open ADAPT_MERGE_V2.pbs
plink --merge-list ADAPT_MergeList.txt --make-bed --out ADAPT_MergedGenotypes_hg19_ATGC
#Error: 16 variants with 3+ alleles present.
# add these 16 variants in the remove list and redo the convert and merge.
cat RemoveTheseAllele.txt ADAPT_MergedGenotypes_hg19_ATGC-merge.missnp >> RemoveTheseAllele_V2.txt
##3.4 Second merge trial 
qsub -A open ADAPT_MERGE_V2.pbs
plink --merge-list ADAPT_MergeList.txt --make-bed --out ADAPT_MergedGenotypes_hg19_ATGC
# since I remove the iid duplicated positions only, there're still rsid have duplicated positions. 
# Remove them in the merge file
grep "Warning" ADAPT_MergedGenotypes_hg19_ATGC.log | awk -F"'" '{print $2,$4}' | tr -s " " "\n" > RemoveTheseDup.txt
rm ADAPT_2784ppl_1M_hg19_ATGC.*
plink --bfile ADAPT_MergedGenotypes_hg19_ATGC --exclude RemoveTheseDup.txt --make-bed --out ADAPT_2784ppl_1M_hg19_ATGC
#1062117 variants loaded from .bim file.
#2784 people (1030 males, 1754 females) loaded from .fam.
#--exclude: 1062021 variants remaining.
#Total genotyping rate is 0.546618.
#1062021 variants and 2784 people pass filters and QC.
#Done for ADAPT

# UIUC2014
## 3.1 Create update list for mismatch positions
# Check if there's any duplicate or mismatch of rsid and pos
cat *.noDI | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk '!a[$0]++' > UIUC2014_ALL.noDI.uniqPOS
awk 'BEGIN{OFS="\t"}{print $1}' UIUC2014_ALL.noDI.uniqPOS | awk '!a[$0]++' | wc -l
#721290
wc -l UIUC2014_ALL.noDI.uniqPOS
#721294 UIUC2014_ALL.noDI.uniqPOS
# 4 duplicate (8) who are they?
# make update list with this
awk -F'\t' 'FNR==NR {x[$1]++; next} x[$1]>1' UIUC2014_ALL.noDI.uniqPOS UIUC2014_ALL.noDI.uniqPOS | head -4 | awk 'BEGIN{OFS="\t"}{print $1,$3}' > UIUC2014_ALL.noDI.uniqPOS.update

##3.2 Create remove list for duplicated iid positions
# check how many have duplicated positions
awk -F'\t' 'FNR==NR {x[$2,$3]++; next} x[$2,$3]>1' UIUC2014_ALL.noDI.uniqPOS UIUC2014_ALL.noDI.uniqPOS > UIUC2014_ALL.noDI.dupPOS
# exctract iid from the dupPOS to get the exclude list
grep "i" UIUC2014_ALL.noDI.dupPOS | awk -F'\t' '{print $1}' > UIUC2014_ALL.noDI.dupPOS.iid
mv UIUC2014_ALL.noDI.dupPOS.iid RemoveTheseAllele.txt
# first merge
chmod u+x UIUC_23aneMeToPlink.bat
./UIUC_23aneMeToPlink.bat
# Error: 70 variants with 3+ alleles present.
# add these 70 variants into remove list, redo the convert and merge.
cat RemoveTheseAllele.txt UIUC2014_MergedGenotypes_hg19_ATGC-merge.missnp >> RemoveTheseAllele.txt
# second merge
chmod u+x UIUC_23aneMeToPlink.bat
./UIUC_23aneMeToPlink.bat
# since I remove the iid duplicated positions only, there're still rsid have duplicated positions. Remove them in the merge file and rename.
grep "Warning" UIUC2014_MergedGenotypes_hg19_ATGC.log | awk -F"'" '{print $2,$4}' | tr -s " " "\n" > RemoveTheseDup.txt
plink --bfile UIUC2014_MergedGenotypes_hg19_ATGC --exclude RemoveTheseDup.txt --make-bed --out UIUC2014_168ppl_705K_hg19_ATGC
#705958 variants loaded from .bim file.
#168 people (76 males, 92 females) loaded from .fam.
#--exclude: 705914 variants remaining.
#Total genotyping rate is 0.814747.
#705914 variants and 168 people pass filters and QC.
# Done for UIUC2014
