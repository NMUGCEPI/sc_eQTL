#####Compute weights, take chr22 as an example
R
wfilename <- '${DIR}/weight/'
bfilename <- '${DIR}/example_snp_chr22'
cfilename <- '${DIR}/example_covar'

options(stringsAsFactors=F)
library(data.table)
library(R.utils)

exp <- read.table('chr22_exp_matrix')
exp <- cbind(ID = rownames(exp), exp)

region <- read.table('chr22_gene_region.xls',h=T)
region <- subset(region, ID %in% colnames(exp))
region$CHR <- region$Chr
region$flank_start <- region$start - 1e6
region$flank_start <- ifelse(region$flank_start<0, 0, region$flank_start)
region$flank_end <- region$end + 1e6

herit <- read.table("Heritability_sigresults_chr22.xls",h=T)
herit <- subset(herit, Vg_p>0.04) 
exp <- exp[,c('ID',herit$gene)]

## start calculation
for(i in 2:ncol(exp)){
	set.seed(1234)
	sub_range <- subset(region, ID==colnames(exp)[i])
	
	write.table(exp[,c(1,1,i)], paste0(wfilename, colnames(exp)[i], ".pheno"), col.names = F, row.names = F, quote = F, sep = "\t")
	commd <- paste0("plink --bfile ", bfilename, " --pheno ", wfilename, colnames(exp)[i], ".pheno --make-bed --out ", wfilename, colnames(exp)[i], " --chr ", sub_range$CHR, " --from-bp ", sub_range$flank_start, " --to-bp ", sub_range$flank_end)
	system(commd)
	rm(list = "commd")
	
	sub_herit <- subset(herit, gene==colnames(exp)[i])
	commd <- paste0("Rscript /data/blj/twas/fusion_twas-master/FUSION.compute_weights.R --bfile ", wfilename, colnames(exp)[i], " --tmp ", wfilename, colnames(exp)[i], ".tmp --covar ", cfilename, " --out ", wfilename, colnames(exp)[i], " --verbose 1 --save_hsq --PATH_gcta /data/blj/twas/fusion_twas-master/gcta_nr_robust --PATH_gemma gemma --PATH_plink plink --models top1,blup,lasso,enet --hsq_set ", sub_herit$Vg_p)
	system(commd)
	rm(list = "commd")
	
	commd <- paste0("rm ", wfilename, colnames(exp)[i], ".b*")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".pheno")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".log")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".fam")
	system(commd)
	rm(list = "commd")
	
	commd <- paste0("rm ", wfilename, colnames(exp)[i], ".nosex")
	system(commd)
	rm(list = "commd")
	
	print(i)
}



