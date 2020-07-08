#!/usr/bin/env Rscript

# for plotting structure of admixture results for each K

#Usage:
#-Place the following files in the same folder as this script
# 	1) your plink binary files (.fam), 
#	2) ADMIXTURE results file (.Q), 
#	3) and 1KG sample info file "i"ntegrated_call_samples_v3.20130502.ALL.panel"

# load file
args = commandArgs(trailingOnly=TRUE)
Qfile<-args[1]  # load .Q file from ADMIXUTRE result
famfile<-args[2] # load fam file from dataset
Popcode1KG<-"integrated_call_samples_v3.20130502.ALL.panel"  # load file for 1KG individual pop info

# get filename and K 
Name <- sub("\\_.*", '', Qfile)
K <- as.numeric(sub('.*?\\.(.*?)\\.Q', '\\1', Qfile))

require(reshape2)
require(ggplot2)

#load the .Q file
qfile=read.table(Qfile,header=F,stringsAsFactors = F)
colnames(qfile)=c(paste0("Anc", 1:K))
#add individual ID labels from .fam file
fam<-read.table(famfile,header=F,stringsAsFactors = F)
colnames(fam)=c("FID","IID","Mot_ID","Fat_ID","Sex","Pheno")
qfile$IID<-fam$IID
# get pop info for 1KG and dataset
Popcode<-read.table(Popcode1KG,header=T)#,sep="\t",stringsAsFactors = F)
qfile$pop<-Popcode$pop[match(qfile$IID,Popcode$sample)]
# change NA into Name
qfile$pop<-as.character(qfile$pop)
qfile$pop[is.na(qfile$pop)]<-Name
#melt the data.frame - i.e. convert to long format
mqfile<-melt(qfile,id.vars=c("IID","pop"))
colnames(mqfile)<-c("IID","pop","Ancestry_component","Ancestry")
#pick one for each super pop for easy read
subpopfile<-mqfile[mqfile$pop %in% c("CHB", "CEU", "YRI","PEL", "GIH", Name),]
# change the order of variable levels with factor()
#mqfile$pop <- factor(mqfile$pop, levels = c("CHB","JPT","CHS","CDX","KHV","CEU","TSI","FIN","GBR","IBS","YRI","LWK","GWD","MSL","ESN","ASW","ACB","MXL","PUR","CLM","PEL","GIH","PJL","BEB","STU","ITU","ADAPT"))

# change the order or IID by anc1
plt<-ggplot(subpopfile,aes(x=factor(IID,levels=unique(IID[order(-Ancestry)])), y=Ancestry, fill=Ancestry_component))+
  geom_bar(stat="identity",width=1,position="fill")

#remove individual IDs from x axis as it can get crowded if you have a lot of samples
plt<-plt+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        panel.background=element_blank())+
  labs(y="Ancestry",fill="Ancestral group") #label axes and legends

#split individuals by population, free_x to show y axis for left column only
plt<-plt+
  facet_wrap(~pop, scale="free_x")#,ncol = 6)

#plt
#one final tweak - change the colors as they are horrendous
#create new pallete to choose from
plt<-plt+
  scale_fill_manual(values=c(
   "skyblue2","#1b9e77","#7570b3","#d95f02", "khaki2","#FB9A99", # lt pink
  "gray70", "yellow3",
  "maroon", "steelblue4",
   "yellow4", "darkturquoise",
  "darkorange4",
  "palegreen2","#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1", "deeppink1", "blue1",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
"dodgerblue2", "#E31A1C", # red
  "green4"))

# save plot as pdf
ggsave(plt,file=paste0(Name,"_K",K,"_admixture.pdf"))


