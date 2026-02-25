##Conduct quality control for scATAC-seq data
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)

library(stringr)
library(rlang)
library(tidyr)
library(purrr)
library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(ggpubr)
library(patchwork)
library(data.table)
library(biovizBase)

scATAC <- readRDS('scATAC.RDS')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
Annotation(scATAC) <- annotations
scATAC <- NucleosomeSignal(scATAC)
scATAC <- TSSEnrichment(scATAC, fast = FALSE)
scATAC1 <- subset(x = scATAC,subset = nCount_ATAC>1000&nCount_ATAC<50000)
scATAC2 <- subset(x = scATAC1,subset = nucleosome_signal < 4)
scATAC <- subset(x = scATAC2,subset = TSS.enrichment > 2)

##Conduct dimensionality reduction and unsupervised clustering
scATAC <- RunTFIDF(scATAC)
scATAC <- FindTopFeatures(scATAC, min.cutoff = 20)
scATAC <- RunSVD(scATAC)
scATAC <- RunUMAP(scATAC, dims = 2:30, reduction = 'lsi',,min.dist = 0.3, n.neighbors = 30L)
scATAC <- FindNeighbors(object = scATAC, reduction = 'lsi', dims = 2:30)
scATAC <- FindClusters(object = scATAC, verbose = FALSE, algorithm = 3)

#Create a gene activity matrix
gene.activities <- GeneActivity(scATAC)
scATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
scATAC <- NormalizeData(
  object = scATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(scATAC$nCount_RNA)
)

##Transfer cell type annotations from the scRNA-seq reference to the scATAC-seq object
scRNA <- readRDS('scRNA.RDS')
DefaultAssay(scATAC) <- 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = scRNA,
  query = scATAC,
  reduction = 'cca')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = scRNA$celltype2,
  weight.reduction = scATAC[['lsi']],
  dims = 2:30)
scATAC@meta.data$celltype <- predicted.labels$predicted.id

##Run MACS2 to call chromatin accessibility peaks
DefaultAssay(scATAC) <- 'ATAC'
peaks <- CallPeaks(
  object = scATAC,
  group.by = "celltype")
# Add cell-type-specific peaks to "ATAC" assay
scATAC@assays$ATAC@misc$celltype_peaks <- peaks

##Visualize gene loci by plotting aggregated scATAC-seq signal
ranges.show <- StringToGRanges(c('chr8-143757708-143757788','chr8-143761931-143762011'))
ranges.show$color <- "#CE0000"
p <- CoveragePlot(
  object = scATAC,
  region = 'PSCA',
  group.by="celltype",
  links = FALSE,
  extend.upstream = 50,
  extend.downstream = 8000,
  ranges = peaks,
  ranges.title = "MACS2",
  region.highlight = ranges.show)
ggsave('CoveragePlot_PSCA_integrated.pdf', p, width = 6, height = 6)


