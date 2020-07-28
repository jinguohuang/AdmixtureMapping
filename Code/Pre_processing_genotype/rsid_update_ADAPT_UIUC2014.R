#!/bin/Rscript

# Pre-processing for genotype data
## update rsid positions with hg19 in ADAPT,UIUC2014 hg18 individuals

# This script can be run by Rscript rsid_update_ADAPT_UIUC2014.R, and will generate 2 update files
# change the file1, bim_file when individual changes

## update iid position in ADAPT, UIUC2014 individuals
file1<-"132001.txt.noDI.i" #change the filename if individual changes
#file1<-"140603.txt.noDI.i" 
#file1<-"140669.txt.noDI.i"
file2<-"UIUC2014_167_noDI_uniqPOS_i"
ONE<-read.table(file1, header=F, sep="\t")
REST<-read.table(file2, header=F, sep="\t")
colnames(ONE)[1:4]<-c("rsid","CHR","POS","AA")
colnames(REST)[1:3]<-c("rsid","CHR","POS")
Merge<-merge(ONE,REST,by=c("rsid"), suffixes = c(".one",".rest"))
length(unique(Merge$rsid))
# are their chrnum same?
Keep<-Merge[as.character(Merge$CHR.one)==as.character(Merge$CHR.rest),]
# update one with new pos
ONE_update<-Keep[c("rsid","CHR.one","POS.rest","AA")]
write.table(ONE_update,paste0(file1,".update"), sep="\t",row.names=F, col.names=F, quote=F)


## update rsid position in ADAPT, UIUC2014 individuals
bim_file<-"132001.txt.noDI.rs" #change the filename if individual changes
#bim_file<-"140603.txt.noDI.rs" 
#bim_file<-"140669.txt.noDI.rs"
ucsc_file<-"UIUC2014_168_rsid_hg19_ucsc"
bim<-read.table(bim_file,sep="\t", header=F)
ucsc<-read.table(ucsc_file,sep="\t", header=F)  
colnames(bim)[1:4] <-c("rsid","CHR","pos","AA")
colnames(ucsc)[1:4] <-c("chrom","pos_start","pos_end","rsid")
# merge 2 table together
Merge<-merge(bim,ucsc,by=c("rsid"))
#compare if pos equals to pos_end
merge_rsid<-unique(Merge$rsid)
# write merge
write.table(Merge,paste0(bim_file,"_Merge"),sep="\t",row.names=F, col.names=F, quote=F)
write.table(merge_rsid,paste0(bim_file,"_Merge_rsid"),sep="\t",row.names=F, col.names=F, quote=F)
Merge$chrom<-sub( "chr","", Merge$chrom)
keep<-Merge[as.character(Merge$CHR)==as.character(Merge$chrom),]
discard<-Merge[!as.character(Merge$CHR)==as.character(Merge$chrom),] 
#reformat the updated position
update<-keep[c("rsid","CHR","pos_end","AA")]
write.table(update,paste0(bim_file,".update"),sep="\t",row.names=F, col.names=F, quote=F)

