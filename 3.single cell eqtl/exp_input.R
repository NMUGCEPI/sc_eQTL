#export the normalized counts for each cell type
library(data.table)
library(dplyr)
library(Seurat)

type <- readRDS("celltype.RDS")
mean_gene_exp <- AverageExpression(type,group.by = 'sample',
slot = 'data') %>% data.frame()
colnames(mean_gene_exp)<-substring(colnames(mean_gene_exp),5)
#Filter the qualified samples of the corresponding cell types
a<-table(type$sample) >= 5
cell_numebr <- as.data.frame(a)
keep_sample <- subset(cell_numebr,a=='TRUE')
gene_counts <- rowSums(GetAssayData(type, assay = "RNA", slot = "counts") > 0)
gene_counts <- as.data.frame(gene_counts)
keep_gene = subset(gene_counts, gene_counts > nrow(type@meta.data)/100)
exp_filter <- mean_gene_exp[rownames(keep_gene),rownames(keep_sample)]
