library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(future)


plan("multiprocess", workers = 12)
#Read in data-------------------------------------
setwd(".")
b2s20 = readRDS("data/merfish_b2s20.rds")
b4s8 = readRDS("data/merfish_b4s8.rds")
b3s20 = readRDS("data/merfish_b3s20.rds")

#Create merged Seurat-------------------------------------
merfish = merge(x = b2s20, y = c(b4s8,b3s20), project = "merfish.merged")

#Integration---------------------------------------------------
# split the dataset into a list of seurat objects (each treatment will be separate object)
options(future.globals.maxSize= 2891289600)
merfish.list <- SplitObject(merfish, split.by = "orig.ident")
# normalize and identify variable features for each library separately
for (i in 1:length(merfish.list)) {
  merfish.list[[i]] <- NormalizeData(merfish.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  #merfish.list[[i]] = FindVariableFeatures(merfish.list[[i]], selection.method = "vst", nfeatures = 3000)
  merfish.list[[i]] <- ScaleData(merfish.list[[i]], features = rownames(merfish.list[[i]]), vars.to.regress = c("nCount_RNA"))
}

# select integration features (all measured genes are used in this case)
merfish.features <- rownames(b2s20)

#Run PCA on each Seurat object
for (name in names(merfish.list)) { 
  merfish.list[[name]] = RunPCA(merfish.list[[name]], features = merfish.features, verbose = TRUE)
}

#Find integration anchors
merfish.anchors <- FindIntegrationAnchors(object.list = merfish.list,
                                                   normalization.method = "LogNormalize", 
                                                   anchor.features = merfish.features,
                                                   reduction = "rpca",
                                                   dims = 1:50) 
# Integrate data, command creates an 'integrated' data assay
all_features <- lapply(merfish.list, row.names)
all_features <- Reduce(intersect, all_features)
merfish <- IntegrateData(anchorset = merfish.anchors,
                                  normalization.method = "LogNormalize",
                                  features.to.integrate = all_features,
                                  dims = 1:50
)


merfish <-  NormalizeData(merfish, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish <- ScaleData(merfish, features = rownames(merfish), assay = "RNA", vars.to.regress = c("nCount_RNA", "orig.ident"))

DefaultAssay(merfish) <- "integrated"
merfish <- ScaleData(merfish, features = rownames(merfish),vars.to.regress = c("nCount_RNA"))
saveRDS(merfish, file = "data/merfish_integrated.rds")

##PCA=====================================================================
merfish <-  RunPCA(merfish, npcs = 100, assay = "integrated")
print(merfish[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish, ndims = 100)
#ggsave(file = "outs/merged/elbow_integrated.pdf", device = "pdf")

VizDimLoadings(merfish, dims = c(1:10), reduction = "pca")
#ggsave(file = "outs/merged/dim_loadings_integrated.pdf", device = "pdf", height = 17, width = 13)

FeatureScatter(merfish, "PC_1","nCount_RNA")
DimPlot(merfish, reduction = "pca", group.by = "orig.ident")
#ggsave(file = "outs/merged/pca_integrated.pdf", device = "pdf")

pdf("outs/merged/Idim_heatmap.pdf", width = 11.69, height = 8.27, useDingbats = F)
DimHeatmap(merfish, dims = c(1:10), cells = 500, balanced = TRUE)
dev.off()
##UMAP/tSNE================================================================
merfish <- RunUMAP(merfish, 
                  dims = 1:30,
                  seed.use = 6,
                  #n.neighbors = 5,
                  #local.connectivity = 15,
                  #assay = "SCT",
                  #n.components = 3
)

DimPlot(merfish, reduction = "umap", pt.size = 0.01, raster = T)
#ggsave(file = "outs/merged/umap_integrated.pdf", device = "pdf", height = 7, width = 10)

##Clustering===============================================================
DefaultAssay(merfish) = "integrated"
merfish = FindNeighbors(merfish,
                       dims = 1:30,
                       reduction = "pca",
                       #k.param = 7,
                       assay = "integrated"
)

merfish = FindClusters(merfish, algorithm = 1, resolution = seq(0.1, 1, by = 0.1))
DimPlot(merfish, reduction = "umap", group.by = c("integrated_snn_res.0.2"), label = T, raster = T, pt.size = 0.1)
#DimPlot(merfish, reduction = "tsne", group.by = c("integrated_snn_res.0.6"), label = T, pt.size = 1)
ggsave(file = "umap_0.2.pdf", device = "pdf", height = 7, width = 9)
#merfish@meta.data[grep(pattern = "^DAM_top100", x = names(merfish@meta.data))] = NULL

clustree(merfish, prefix ="integrated_snn_res.")
#ggsave(filename = "outs/merged/clustree.pdf", device = "pdf", height = 20)

FeaturePlot(merfish, c("nCount_RNA", "nFeature_RNA","volume", "blank_sum", "blank_mean"),
                 reduction = "umap", pt.size = 0.1, min.cutoff = "q1", max.cutoff = "q99")
#ggsave(file = "outs/merged/QC.pdf", device = "pdf", height = 10, width = 10)

VlnPlot(merfish, c("nCount_RNA", "nFeature_RNA","volume", "blank_sum", "blank_mean"),
               pt.size = 0, group.by = "integrated_snn_res.0.2")
#ggsave(file = "outs/merged/QC_vln_0.2.pdf", device = "pdf", height = 8, width = 14)

#Read in marker gene sets
consmarkerList = readRDS("~/GOKCE LAB/Code/gene_sets/consensusMarkerList.rds")

#Plot known marker gene set expression
cellmark_merfish = list()
for (celltype in names(consmarkerList)){
  cellmark_merfish[[celltype]] = intersect(rownames(merfish), consmarkerList[[celltype]])
}

cellmark_merfish = cellmark_merfish[lapply(cellmark_merfish,length)>0]
for (celltype in names(cellmark_merfish)){
  FeaturePlot(merfish, features = cellmark_merfish[[celltype]], ncol = 3, reduction = "umap", pt.size = 0.1)
  ggsave(file = paste0( celltype, ".pdf"), device = "pdf", height = 16, width = 10)
}

# Find markers==============================================================
Idents(merfish) <- "integrated_snn_res.0.2"
DefaultAssay(merfish) = "RNA"
markers_0.2 = FindAllMarkers(merfish,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             
write.table(markers_0.2, file = "markers_0.2.txt", quote = F, sep = "\t", col.names = NA)

# Plot cluster markers
to.plot = markers_0.2 %>% group_by(cluster) %>% top_n(5,avg_log2FC) %>% as.data.frame()
cols = rev(colorRampPalette(brewer.pal("RdBu",n = 11))(256))
for(i in levels(unique(Idents(merfish)))){
  FeaturePlot(merfish, features = as.character(to.plot[to.plot$cluster == i, "gene"]), pt.size = 0.3, ncol = 4)
  ggsave(paste0("top5markers_cluster",i,".png"), device = "png", width = 13, height = 10)
}

saveRDS(merfish, file = "data/merfish_integrated.rds")




































