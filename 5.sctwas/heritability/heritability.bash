##GCTA calculate GRM, take chr22 as an example
for gene in $(cat genelist)
do
plink --noweb --bfile example_snp_chr22 --extract gene_region/$gene --range --make-bed --out gene_plink/$gene
/data/gc_sceqtl/herit/GCTA/gcta64  --bfile gene_plink/$gene  --autosome  --make-grm  --out gene_grm/$gene
/data/gc_sceqtl/herit/GCTA/gcta64  --reml --grm gene_grm/$gene --pheno pheno/$gene --covar ./ccovar --qcovar ./qcovar --reml-no-constrain --out heritability/$gene
done

##Summarize the heritability data
cd ./heritability
ls -1 | grep ".hsq" > content

R
gene<-read.table('./content',stringsAsFactors=F)
n<-nrow(gene)
heri<-data.frame(gene=1:n,Vg=1:n,SEg=1:n,Ve=1:n,SEe=1:n,Vp=1:n,SEp=1:n,Vg_p=1:n,SEg_p=1:n,
		logL=1:n,logLref=1:n,LRT=1:n,df=1:n,P=1:n,Num=1:n)
for(i in 1:n){
	data<-read.table(gene[i,1],h=T,stringsAsFactors=F,fill=T)
	heri$gene[i]<-gene[i,1]
	heri$Vg[i]<-data[1,2]
	heri$SEg[i]<-data[1,3]
	heri$Ve[i]<-data[2,2]
	heri$SEe[i]<-data[2,3]
	heri$Vp[i]<-data[3,2]
	heri$SEp[i]<-data[3,3]
	heri$Vg_p[i]<-data[4,2]
	heri$SEg_p[i]<-data[4,3]
	heri$logL[i]<-data[5,2]
	heri$logLref[i]<-data[6,2]
	heri$LRT[i]<-data[7,2]
	heri$df[i]<-data[8,2]
	heri$P[i]<-data[9,2]
	heri$Num[i]<-data[10,2]
}
heri$gene<-substr(heri$gene,1,15)
heri_sig<-subset(heri,P<0.05 & Vg>0 & Vp>0)
write.table(heri_sig,file='Heritability_sigresults_chr22.xls',row.names=F,quote=F)
