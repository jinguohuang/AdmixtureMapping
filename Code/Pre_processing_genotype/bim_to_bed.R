#!/bin/Rscript

# Pre-processing for genotype data - CV

# convert bim file to bed that can use LiftOver
# This script need input bim file, will generate a bed file

BimToBed<-function(bim){
  #load file
  bim_file<-read.table(bim,sep="\t",header=F, colClasses = c("numeric" , "character" , "numeric" , "numeric" , "character" , "character"),comment.char = "", check.names = FALSE)
  colnames(bim_file)[1:6]<-c("CHR","rsid","cm","pos","A1","A2")
  # add new colomn "chrom"
  bim_file$chrom<-bim_file$CHR
  # change 23 to X
  bim_file$chrom<-gsub("23","X",bim_file$chrom)
  # change 24 to Y
  bim_file$chrom<-gsub("24","Y",bim_file$chrom)
  # change 25 to XY
  bim_file$chrom<-gsub("25","X",bim_file$chrom)
  # change 26 to MT
  bim_file$chrom<-gsub("26","M",bim_file$chrom)
  # add "chr" to the chrom
  bim_file$chrom<-sub("^", "chr", bim_file$chrom)
  #disabling scientific notation in R, avoid big number show in expanational way
  options(scipen = 999)
  # add new column "chromStart" the same with pos-1
  bim_file$chromStart<-bim_file$pos-1
  # add new column "chromEnd" the same with pos
  bim_file$chromEnd<-bim_file$pos
  # re order by column name
  bed_file <- bim_file[c("chrom", "chromStart", "chromEnd", "rsid")]

  # write.table
  write.table(bed_file, paste0(bim,"_bimtobed"), row.names=F, col.names=F, quote=F, sep="\t")
}

args <- commandArgs(trailingOnly = TRUE)
bim<-args[1]
BimToBed(bim)
