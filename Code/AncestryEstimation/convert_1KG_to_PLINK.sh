#!/bin/bash

# Ancestry estimation
## Convert 1KG into Plink1 format

# download 1KG plink file https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
# download 3 bolded files as suggested on PLINK website
# Merged datasets,excludes extra chrM samples
wget https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
wget https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1
# download common sample info
wget https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1

# rename so PLINK2 can recognize files
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
mv 'phase3_corrected.psam?dl=1' all_phase3.psam

# install PLINK2 to decompress zst files
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20200607.zip
unzip plink2_linux_avx2_20200607.zip
# make plink2 executable
chmod +x plink2

# decompress zst file 
./plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
# Then convert PLINK2 binary format to PLINK1 binary format
./plink2 --pfile all_phase3 vzs --max-alleles 2 --make-bed --out all_phase3


