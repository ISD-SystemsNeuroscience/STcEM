library(Seurat)
library(data.table)
library(stringr)
library(Matrix)
library(data.table)
library(mltools)
library(tidyverse)
library(RColorBrewer)
library(clipr)

setwd(".")

#Read Allen Brain Atlas reference data
expression_ABA = fread("/Users/pandrovi/GOKCE LAB/Code/reference_datasets/Allen_brain_20211220/data/matrix.csv", data.table = FALSE) #large file, use data.table::fread() to read it faster
metadata = fread("/Users/pandrovi/GOKCE LAB/Code/reference_datasets/Allen_brain_20211220/data/metadata.csv", data.table = FALSE)
expression_ABA = Matrix(as.matrix(expression_ABA), sparse = TRUE)
str(metadata)
expression_ABA[1:10,1:10]

#Split it into major cell classes (Glut, GABA, non-neuronal (not used))
## Downsample GABAergic further to avoid memory limits, full dataset is prohibitively large
meta_GABA = filter(metadata, class_label == "GABAergic")
plot(table(meta_GABA$subclass_label))
#Downsample glutamatergic cells based on label to cap the max number of cells per subclass
meta_GABA_downsampled <- meta_GABA %>% group_by(subclass_label) %>% slice_sample(n=15000) #max 15000 cells per subclass
plot(table(meta_GABA_downsampled$subclass_label))
expression_GABA = filter(expression_ABA, sample_name %in% meta_GABA_downsampled$sample_name)
rownames(expression_GABA) = expression_GABA$sample_name
expression_GABA$sample_name = NULL
expression_GABA[1:10,1:10]
expression_GABA = Matrix(as.matrix(expression_GABA), sparse = TRUE)

## Downsample Glutamatergic further to avoid memory limits, full dataset is prohibitively large
meta_glut = filter(metadata, class_label == "Glutamatergic")
plot(table(meta_glut$subclass_label))
#Downsample glutamatergic cells based on label to cap the max number of cells per subclass
meta_downsampled <- meta_glut %>% group_by(subclass_label) %>% slice_sample(n=15000) #max 15000 cells per subclass
plot(table(meta_downsampled$subclass_label))

table(meta_downsampled$subclass_label) - table(meta_glut$subclass_label) #Which subclassess were downsampled?
expression_Glut = filter(expression_ABA, sample_name %in% meta_downsampled$sample_name)
rownames(expression_Glut) = expression_Glut$sample_name
expression_Glut$sample_name = NULL
expression_Glut[1:10,1:10]
expression_Glut = Matrix(as.matrix(expression_Glut), sparse = TRUE)

#saveRDS(expression_non_neuronal, file = "non_neuronal.rds")
saveRDS(expression_GABA, file = "data/ABA_GABAergic_downsampled.rds")
saveRDS(expression_Glut, file = "data/ABA_Glutamatergic_downsampled.rds")

#Create Seurat object, add metadata
expression_GABA = t(expression_GABA)
expression_Glut = t(expression_Glut)
table(rownames(expression_GABA) == rownames(expression_Glut))
expression_neuronal = cbind(expression_GABA, expression_Glut)
ABA_neuronal = CreateSeuratObject(counts = expression_neuronal, project = "ABA_neuronal", min.cells = 20, min.features = 200)
metadata = metadata %>% filter(sample_name %in% colnames(ABA_neuronal))
rownames(metadata) = metadata$sample_name
ABA_neuronal = AddMetaData(ABA_neuronal, metadata = metadata)
table(ABA_neuronal$subclass_label)

#saveRDS(ABA_neuronal, "data/Seurat_ABA_Neuronal_Downsampled.rds")

# Data normalization, scaling and variable feature selection=======
ABA_neuronal[["percent.mt"]] <- PercentageFeatureSet(ABA_neuronal, pattern = "^mt-")
VlnPlot(ABA_neuronal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
VlnPlot(ABA_neuronal, group.by = "subclass_label", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ABA_neuronal <- subset(ABA_neuronal, percent.mt <= 15)

ABA_neuronal <-  NormalizeData(ABA_neuronal, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
ABA_neuronal <- FindVariableFeatures(ABA_neuronal, selection.method = "vst", nfeatures = 3000)
ABA_neuronal <- ScaleData(ABA_neuronal, features = VariableFeatures(ABA_neuronal), vars.to.regress = "nCount_RNA")

top30 <- head(VariableFeatures(ABA_neuronal), 30)
plot1 <- VariableFeaturePlot(ABA_neuronal, pt.size = 0.1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = T)
plot2
ggsave(plot2, file = "outs/var_features.pdf", device = "pdf")


#Check correlation of sequencing depth and expression of some genes in original and normalized data
FeatureScatter(ABA_neuronal, "nCount_RNA","Map2")
FeatureScatter(ABA_neuronal, "nCount_RNA", "Cd74")
FeatureScatter(zalohaABA_neuronal, "nCount_RNA","Cx3cr1")
FeatureScatter(zalohaABA_neuronal, "nCount_RNA", "Cd74")

##PCA=====================================================================
ABA_neuronal <-  RunPCA(ABA_neuronal, npcs = 100)
print(ABA_neuronal[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(ABA_neuronal, ndims = 100)
ggsave(file = "outs/elbow.pdf", device = "pdf")

VizDimLoadings(ABA_neuronal, dims = c(1:10), reduction = "pca")
ggsave(file = "outs/dim_loadings.pdf", device = "pdf", height = 17, width = 11)

FeatureScatter(ABA_neuronal, "PC_1","nCount_RNA")
DimPlot(ABA_neuronal, reduction = "pca", group.by = "Group")
ggsave(file = "outs/pca.pdf", device = "pdf")

pdf("outs/dim_heatmap.pdf", width = 11.69, height = 8.27, useDingbats = F)
DimHeatmap(ABA_neuronal, dims = c(1:10), cells = 500, balanced = TRUE)
dev.off()

##UMAP/tSNE================================================================
ABA_neuronal <- RunUMAP(ABA_neuronal, 
                  dims = 1:30,
                  seed.use = 6,
                  #n.neighbors = 5,
                  #local.connectivity = 15,
                  #assay = "SCT",
                  #n.components = 3
)
ABA_neuronal <- RunTSNE(ABA_neuronal, dims = 1:10) 
#colors = c("#b2182b","#fdae61","#2166ac","#92c5de")
cols_subclass = unique(ABA_neuronal$subclass_color)
names(cols_subclass) = unique(ABA_neuronal$subclass_label)
DimPlot(ABA_neuronal, reduction = "umap", group.by = "subclass_label", cols = cols_subclass, pt.size = 0.1, raster = FALSE)
ggsave(file = "outs/umap_30dims.pdf", device = "pdf", height = 7, width = 10)

DimPlot(ABA_neuronal, reduction = "umap", group.by = "cortical_layer_order", pt.size = 0.1)
#ggsave(file = "outs/tsne_groups_10dims.pdf", device = "pdf", height = 7, width = 10)

# Add umap provided by Allen Brain Institute (optional)
umap_ABA = read.csv(file = "data/tsne.csv")
head(umap_ABA)
colnames(umap_ABA) = c("sample_name", "umap_ABA_1", "umap_ABA_2")
cells.to.keep = intersect(umap_ABA$sample_name, colnames(ABA_neuronal))
umap_ABA = umap_ABA %>% filter(sample_name %in% cells.to.keep)
rownames(umap_ABA) = umap_ABA$sample_name
umap_ABA$sample_name = NULL
umap_ABA_assay = CreateDimReducObject(embeddings = as.matrix(umap_ABA), key = "umapABA_")
ABA_neuronal@reductions[["umapABA"]] = umap_ABA_assay

DimPlot(ABA_neuronal, reduction = "umapABA", group.by = "subclass_label", cols = cols_subclass, pt.size = 0.1, raster = FALSE)
#ggsave(file = "outs/umap_ABA.pdf", device = "pdf", height = 7, width = 10)

#Save Seurat object
saveRDS(ABA_neuronal, "data/Seurat_ABA_Neuronal_Downsampled.rds")

