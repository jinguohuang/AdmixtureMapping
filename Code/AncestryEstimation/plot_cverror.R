#!/usr/bin/env Rscript
# for plotting crossvalidation error of admixture results

#Usage:
#	put cverror file in the same folder with this script 

# load file
args = commandArgs(trailingOnly=TRUE)
cverror<-args[1]

# get filename 
Name <- sub("\\_.*", '', cverror)
library(ggplot2)
cverror<-read.table(cverror)
colnames(cverror)<-c("K","cv_error")

# plot cv_error against K
plt<-ggplot(data=cverror, aes(x=K, y=cv_error, group=1)) + 
geom_line() + geom_point()+ 
xlab("K") + ylab("Cross-validation error") + 
ggtitle(paste0("Cross validation error in ", Name))
ggsave(plt,file=paste0(Name,"_cverror.pdf"))

