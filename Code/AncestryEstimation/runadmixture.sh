#!/bin/bash

# Ancestry estimation
## Ancestry estimation using ADMIXTURE

# install ADMIXTURE
thisdir=$(pwd)
wget http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xvzf admixture_linux-1.3.0.tar.gz

# get filename
ls *.bed |cut -d'.' -f1 > ALL

for file in $(cat ALL)
do
	for k in {2..16}
	do
		echo "#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=20:00:00
#PBS -l pmem=32gb
#PBS -A tll30_a_g_bc_default 
#PBS -j oe
cd ${thisdir}
admixture --cv -j6 ${file} ${k}" >> job_${file}_k${k}.pbs #Creating job to run admixture
	qsub job_${file}_k${k}.pbs # submit job to ICS to run
	done
done

