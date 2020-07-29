#!/bin/Rscript

# Pre-processing for genotype data - CV
# This script will generate a ucsc hg19 position file of CV, an update position map, and a remove list

##0 Get hg19 position for all variants
# upload rsid list to UCSC genome browser to get hg19 positions
# due to the large number of rsid, if upload to dbSNP151 all, the website will return error
# so upload to dbSNP151 common first, upload those who cannot find match in common to db151all
# Combine these files to make the hg19 position file
common<-read.table("CV_151common", header=F, sep="\t") #929984
aa<-read.table("CV_151all_aa", header=F, sep="\t") #7733
ab1<-read.table("CV_151all_ab_aa", header=F, sep="\t") #3859
ab2<-read.table("CV_151all_ab_ab_2aa", header=F, sep="\t") #1909
ab3<-read.table("CV_151all_ab_ab_2ab", header=F, sep="\t") #9288
ac<-read.table("CV_151all_ac", header=F, sep="\t") #7531
ad<-read.table("CV_151all_ad", header=F, sep="\t") #7599
ae<-read.table("CV_151all_ae", header=F, sep="\t") #7748
af<-read.table("CV_151all_af", header=F, sep="\t") #3410
all<-rbind(common,aa,ab1,ab2,ab3,ac,ad,ae,af) #929984
#remove duplicate
all_1<-all[!grepl("_", all$V1),] #905243
dup<-all_1[duplicated(all_1$V4)|duplicated(all_1$V4,fromLast = TRUE),] #1752, chrX and chrY, position are different
write.table(all_1, "CV_rsid_hg19_ucsc",sep="\t",row.names=F, col.names=F, quote=F)

#1 Update positions that have same chromosome number
# merge and update those with same chr num
# load file
bim_file<-"CV_hg18.bim.rs"
ucsc_file<-"CV_rsid_hg19_ucsc"
bim<-read.table(bim_file,sep="\t", header=F) #906322
ucsc<-read.table(ucsc_file,sep="\t", header=F) #905243 
colnames(bim)[1:6] <-c("CHR","rsid","cm","pos","A1","A2")
colnames(ucsc)[1:4] <-c("chrom","pos_start","pos_end","rsid")
# merge 2 table together
Merge<-merge(bim,ucsc,by=c("rsid")) #905243
#compare if pos equals to pos_end
merge_rsid<-unique(Merge$rsid) #904367
dup1<-Merge[duplicated(Merge$rsid)|duplicated(Merge$rsid,fromLast = TRUE),]
#change the CHR into ucsc format
Merge$CHR<-sub( "23","X", Merge$CHR)
Merge$CHR<-sub( "24","Y", Merge$CHR)
Merge$CHR<-sub( "25","X", Merge$CHR)
Merge$CHR<-sub( "^","chr", Merge$CHR)
#keep those chromosome number same
keep<-Merge[as.character(Merge$CHR)==as.character(Merge$chrom),] #904312
keep_rsid<-unique(keep$rsid) #904312
write.table(keep,paste0(bim_file,".merge_keep"),sep="\t",row.names=F, col.names=F, quote=F)
# make update list: rsid and position
update<-keep[c("rsid","pos_end")]
write.table(update,paste0(bim_file,".update"),sep="\t",row.names=F, col.names=F, quote=F)
# who's not match?
discard<-Merge[!as.character(Merge$CHR)==as.character(Merge$chrom),] #931
# except chrX and chrY, whose chr num changed?
discard_auto<-discard[!grepl("chrX|chrY", discard$CHR),] #50
#2 remove list come up with those who are in keep not in bim
remove<-bim[!bim$rsid %in% keep$rsid,] #2010
write.table(remove,paste0(bim_file,".removed"),sep="\t",row.names=F, col.names=F, quote=F)