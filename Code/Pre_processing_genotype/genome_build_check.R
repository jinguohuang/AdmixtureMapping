#!/bin/Rscript

# Pre-processing for genotype data
## Genome build check

# Prepare 2 files for this script： bim file and ucsc position file
# This script will come back to 3 files with matched and unmatched record

# create function for the genome build check:
GenomeBuild<-function(bim,ucsc){  
  # load files
  bim<-read.table(bim_file,sep="\t", header=F,comment.char = "", check.names = FALSE) # '#' not comment out in this file
  ucsc<-read.table(ucsc_file,sep="\t", header=F)
  # rename files   
  colnames(bim)[1:6] <-c("CHR","rsid","BP","pos","A1","A2")
#  colnames(bim)[1:4] <-c("CHR","rsid","BP","pos") # rename if load with 4 column bim file
  colnames(ucsc)[1:4] <-c("chrom","pos_start","pos_end","rsid")
  # merge 2 table together
  Merge<-merge(bim,ucsc,by=c("rsid"))
  #compare if pos equals to pos_end
  merge_rsid<-unique(Merge$rsid)
  # get matched position where pos==pos_end
  match<-Merge[Merge$pos==Merge$pos_end, ]
  # get matched rsid no duplication
  match_rsid<-unique(match$rsid)
  # get unmatched position
  unmatch<-Merge[!Merge$rsid %in% match_rsid, ]
  unmatch_rsid<-unique(unmatch$rsid)
  # write merged file
  write.table(Merge,paste0(bim_file,"_Merge"),sep="\t",row.names=F, col.names=T, quote=F)
  write.table(merge_rsid,paste0(bim_file,"_Merge_rsid"),sep="\t",row.names=F, col.names=F, quote=F)
  # write matched position
  write.table(match,paste0(bim_file,"_Merge_match"),sep="\t",row.names=F, col.names=T, quote=F)
  write.table(match_rsid,paste0(bim_file,"_Merge_match_rsid"),sep="\t",row.names=F, col.names=F, quote=F)
  # write unmatched position
  write.table(unmatch,paste0(bim_file,"_Merge_unmatch"),sep="\t",row.names=F, col.names=T, quote=F)
  write.table(unmatch_rsid,paste0(bim_file,"_Merge_unmatch_rsid"),sep="\t",row.names=F, col.names=F, quote=F)
}

#bim_file <- "Euro180_176ppl_312K_hg19.bim"
#ucsc_file <- "Euro180_rsid_db153"
args <- commandArgs(trailingOnly = TRUE)
# input files
bim_file<-args[1] 
ucsc_file<-args[2]
# call function to get the output files
GenomeBuild(bim,ucsc)