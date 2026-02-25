##Conduct single cell quality control, dimensionality reduction and unsupervised clustering.
library(Seurat)
data <- Read10X(data.dir = 'input') #The variable "input" represents the directory path where the matrix files are stored.
data <- CreateSeuratObject(counts = data ,project = 'seurat', min.cells = 3, min.features = 200,names.delim = '_')
data[["percent.mt"]]<-PercentageFeatureSet(object = data, pattern = "^MT-")
data<-subset(x=data,subset=nFeature_RNA>200&nFeature_RNA<4000&percent.mt<20&nCount_RNA>1000&nCount_RNA<20000) 
data <- NormalizeData(object = data, normalization.method = 'LogNormalize',scale.factor = 10000)
data <- FindVariableFeatures(object = data,selection.method='vst',mean.function = ExpMean,dispersion.function=LogVMR,mean.cutoff=c(0.125,3),dispersion.cutoff=c(0.5,Inf))
data <- ScaleData(data,vars.to.regress = c("percent.mt","nCount_RNA","orig.ident"))
data <- RunPCA(object=data,npcs=100,pc.genes=VariableFeatures(object=data))
data <- FindNeighbors(object=data,dims=1:50)
data <- FindClusters(object=data,resolution=0.3)
data <- RunUMAP(data,dims = 1:50,min.dist = 0.3, n.neighbors = 25L)
