##!/bin/bash

# Ancestry Estimation
## Plot for ADMIXTURE results

# make filename list
$ ls *.fam | cut -d'.' -f1,2 > ALL

# Plot population structure for each dataset and each K
for file in $(cat ALL)
do
	for k in {2..16}
	do
		echo "plot for $file with K=$k"
		Rscript plot_admixture.R ${file}.${k}.Q ${file}.fam
	done
done

# Create file that extract crossvalidation error and plot it
for name in $(cat ALL)
do
	echo "Extracting crossvalidation error for $name"
	grep -h CV job1_${name}*.pbs.o*| cut -d'=' -f2 | awk -F')' '{print $1,$2}'| awk -F":" '{print $1,$2}' | sort -k1n > ${name}_cverror
	echo "Plotting crossvalidation error for $name" 
	Rscript plot_cverror.R ${name}_cverror	
done

