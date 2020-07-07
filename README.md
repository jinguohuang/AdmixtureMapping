# Admixture Mapping Project

Here is the pipeline for the admxiture mapping analysis.

## Ancestry Estimation Pipeline
1) Strain check with snpflip
2) Harmonize datasets with reference data using Genotype Harmonizer
3) Convert 1KG to PLINK1 file format
4) Merge our datasets with 1KG and QC
5) LD based SNP pruning
6) Ancestry estimation using ADMIXTURE
7) Plot for ADMIXTURE results
