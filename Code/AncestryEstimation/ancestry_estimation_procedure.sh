#!/bin/bash

# Ancestry Estimation
## Procedures to run each scripts 

# 1) Strain check with snpflip
chmod +x strain_check_snpflip.sh
./strain_check_snpflip.sh

# 2) Harmonize datasets with reference data using Genotype Harmonizer
chmod +x harmonize_with_1KG.sh
./harmonize_with_1KG.sh

# 3) Convert 1KG to PLINK1 file format
chmod +x convert_1KG_to_PLINK.sh
./convert_1KG_to_PLINK.sh

# 4) Merge our datasets with 1KG and QC
chmod +x merge_with_1KG.sh
./merge_with_1KG.sh

# 5) LD based SNP pruning
chmod +x QC_ld_pruning.sh
./QC_ld_pruning.sh

# 6) Ancestry estimation using ADMIXTURE
chmod +x runadmixture.sh
./runadmixture.sh

# 7) Plot for ADMIXTURE results
chmod +x plot_admixture_cverror.sh
./plot_admixture_cverror.sh


