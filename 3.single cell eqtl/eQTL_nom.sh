#!/bin/bash

#Use FastQTL to analysis eqtl
run_FastQTL_threaded.py \
sample.vcf.gz\
celltype_exp.bed.gz\
celltype_eqtl\
--covariates covariates.txt \
--window 1e6 --chunks 100 --threads 16 \
--output_dir output_path

#Use tensorqtl to analysis eqtl
pip3 install tensorqtl
#genetype output
plink2 \
    --vcf sample.vcf.gz \
    --out sample
#eqtl analysis	
python3 -m tensorqtl \
    ${GenoDir}/sample \
    ${DataDir}/celltype_exp.bed.gz \
    ${OUT_DIR} \
    --covariates covariates.txt \
    --mode cis_nominal
