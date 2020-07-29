#!/bin/Rscript

# Pre-processing for genotype data - Euro180

##0 Get hg19 position for all variants
# upload rsid list to UCSC genome browser to get hg19 positions
# due to the large number of rsid, if upload to dbSNP151 all, the website will return error
# so upload to dbSNP151 common first, upload those who cannot find match in common to db151all
# Combine these 2 files to make the hg19 position file
file1<-"Euro180_db151common"
file2<-"Euro180_db151all"
common<-read.table(file1,header = F,sep="\t")
all<-read.table(file2,header = F,sep="\t")
db151<-rbind(common,all)
# delete duplicate
colnames(db151)[1:4] <-c("chrom","pos_start","pos_end","rsid")
a<-db151$rsid #320783
# ALL the dupplicated rows
b<-db151[duplicated(a)|duplicated(a,fromLast = TRUE),]
# who is duplicated
c<-unique(b$rsid)
# delete chromosome with "_" 
d<-db151[!grepl("_", db151$chrom),]#314391
# save as file
write.table(d,"Euro180_ucsc_hg19",sep="\t",row.names=F, col.names=F, quote=F)
length(unique(d$rsid))
#314389
e<-d[duplicated(d$rsid)|duplicated(d$rsid,fromLast = TRUE),]
#4 duplicate are 2 rsid have both chrX and chrY

##1. Prepare position update list
# Get update list. Merge the bim with ucsc file, keep those have same chromosome number. Discard the rest of them.
# load the bim file and the ucsc position file
bim<-read.table(bim_file,sep="\t", header=F) # ,comment.char = "", check.names = FALSE) # '#' not comment out in this file
ucsc<-read.table(ucsc_file,sep="\t", header=F)
# change colnames  
colnames(bim)[1:6] <-c("CHR","rsid","BP","pos","A1","A2")
colnames(ucsc)[1:4] <-c("chrom","pos_start","pos_end","rsid")
# merge 2 table together
Merge<-merge(bim,ucsc,by=c("rsid"))
# write merge file
write.table(Merge,paste0(bim_file,"_Merge"),sep="\t",row.names=F, col.names=F, quote=F)
# keep those CHR==chrom
file<-"Euro180_176ppl_317K_hg19_ATGC.bim_Merge"
Merge<-read.table(file, header=F, sep="\t")
colnames(Merge)[1:9]<-c("rsid","CHR","BP","pos","A1","A2","chrom","pos_start","pos_end")
# reformat chrom into CHR, bim file and ucsc have different chromosome format
Merge$chrom<-sub( "chr","", Merge$chrom)
Merge$chrom<-sub( "X","23", Merge$chrom)
Merge$chrom<-sub( "Y","24", Merge$chrom)
# keep those positions that bim and ucsc are on the same chromosome, discard those not
keep<-Merge[as.character(Merge$CHR)==as.character(Merge$chrom),] #314383
discard<-Merge[!as.character(Merge$CHR)==as.character(Merge$chrom),] #8
# create the update rsid map for Plink --update-map
update<-keep[c("rsid","pos_end")] #314383
write.table(update,"Euro180_update_map",sep="\t",row.names=F, col.names=F, quote=F)

##2 Prepare exclude list
# Get remove list by exlude all update positions from bim. 
bim<-"Euro180_176ppl_317K_hg19_ATGC.bim"
bim_file<-read.table(bim, header=F, sep="\t") #317503
remove<- bim_file[!bim_file$V2 %in% update$rsid,] #3120
remove_rsid<-remove$V2
write.table(remove_rsid,"Euro180_remove_rsid",sep="\t",row.names=F, col.names=F, quote=F)