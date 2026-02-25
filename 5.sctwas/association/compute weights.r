###compute association
##Use fusion tools http://gusevlab.org/projects/fusion/
###Download and unpack the FUSION software package from github:
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
cd fusion_twas-master
###Download and unpack the (1000 Genomes)  LD reference data:
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2

Rscript ${DIR}/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats ${DIR}/example_gwas.sumstats \
    --weights ${DIR}/example_weights_file.list \
    --weights_dir ${DIR}/weight \
    --ref_ld_chr ${DIR}/LDREF_EAS/1000G.EAS. \
    --chr 22 \
    --GWASN 21168 \
    --min_r2pred 0.4 \
    --out ${DIR}/twas.chr22.dat
done