
setwd(".")

#Read in data
merfish = readRDS("data/merfish_integrated.rds")
ABA_neuronal = readRDS("data/Seurat_ABA_Neuronal_Downsampled.rds")


# NEURO ONLY -  Subset integrated merfish data to neuronal cells only, and re-run--------------------------------
neuronal_clusters = c(0,4,8,10, 14) #Collect numbers of neuronal clusters from resolution 0.2. Exclude striatal spiny neuron cluster 2 (not included in cortex+hippocampus reference)
merfish_neuro = subset(merfish, integrated_snn_res.0.2 %in% neuronal_clusters)
DimPlot(merfish_neuro, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

merfish_neuro <-  NormalizeData(merfish_neuro, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_neuro <- ScaleData(merfish_neuro, features = rownames(merfish_neuro), vars.to.regress = "nCount_RNA")

top30 <- head(VariableFeatures(merfish_neuro), 30)
plot1 <- VariableFeaturePlot(merfish_neuro, pt.size = 0.1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = T)
plot2
#ggsave(plot2, file = "outs/var_features.pdf", device = "pdf")

##PCA=====================================================================
merfish_neuro <-  RunPCA(merfish_neuro, npcs = 100, assay = "integrated")
print(merfish_neuro[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish_neuro, ndims = 100)
ggsave(file = "outs/merged/elbow_integrated_neuro.pdf", device = "pdf")

VizDimLoadings(merfish_neuro, dims = c(1:10), reduction = "pca")
ggsave(file = "outs/merged/dim_loadings_integrated_neuroOnly.pdf", device = "pdf", height = 17, width = 13)

FeatureScatter(merfish_neuro, "PC_1","nCount_RNA")
DimPlot(merfish_neuro, reduction = "pca", group.by = "orig.ident")
ggsave(file = "outs/merged/pca_integrated_neuroOnly.pdf", device = "pdf")

pdf("outs/merged/Idim_heatmap.pdf", width = 11.69, height = 8.27, useDingbats = F)
DimHeatmap(merfish_neuro, dims = c(1:10), cells = 500, balanced = TRUE)
dev.off()

##UMAP/tSNE================================================================
merfish_neuro <- RunUMAP(merfish_neuro, 
                         dims = 1:25,
                         seed.use = 6,
                         #n.neighbors = 5,
                         #local.connectivity = 15,
                         #assay = "SCT",
                         #n.components = 3
)
merfish_neuro <- RunTSNE(merfish_neuro, dims = 1:10) 
#colors = c("#b2182b","#fdae61","#2166ac","#92c5de")
cols_subclass = unique(merfish_neuro$subclass_color)
names(cols_subclass) = unique(merfish_neuro$subclass_label)
DimPlot(merfish_neuro, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, raster = FALSE)
ggsave(file = "outs/umap_25dims_neuroOnly.pdf", device = "pdf", height = 7, width = 10)

DimPlot(merfish_neuro, reduction = "umap", group.by = "cortical_layer_order", pt.size = 0.1)
ggsave(file = "outs/tsne_groups_10dims.pdf", device = "pdf", height = 7, width = 10)

##Transfer labels--------------------------------------------------
options(future.globals.maxSize= 2891289600)
common.genes = intersect(rownames(merfish_neuro), rownames(ABA_neuronal))
common.genes
merfish_neuro_anchors = FindTransferAnchors(reference = ABA_neuronal, query = merfish_neuro, normalization.method = "LogNormalize", reduction = "rpca",
                                            #reference.assay = "integrated",  #reference.reduction = "pca",
                                            features = common.genes, dims = 1:25, k.anchor = 5)
#saveRDS(object = merfish_neuro_anchors, file = "merfish_neuro_neuronalABA.anchors.rds")

predictions.neuro = TransferData(anchorset = merfish_neuro_anchors, refdata = ABA_neuronal$class_label,
                                 #weight.reduction = "cca",
                                 dims = 1:25)
merfish_neuro@meta.data[colnames(predictions.neuro)] = NULL
head(predictions.neuro)
colnames(predictions.neuro)[c(1,4)] = c("class_label", "prediction.score.max.class_label")
merfish_neuro= AddMetaData(merfish_neuro, metadata = predictions.neuro)

predictions.subclass = TransferData(anchorset = merfish_neuro_anchors, refdata = ABA_neuronal$subclass_label,
                                    #weight.reduction = "cca",
                                    dims = 1:25)
head(predictions.subclass)
merfish_neuro@meta.data[colnames(predictions.subclass)] = NULL
colnames(predictions.subclass)[c(1,38)] = c("subclass_label", "prediction.score.max_subclass")
merfish_neuro= AddMetaData(merfish_neuro, metadata = predictions.subclass)

head(merfish_neuro@meta.data)


saveRDS(predictions.neuro, file = "predictions_fromNeuroToNeuro.rds")

## Visualize results----------------------------------------------
DimPlot(merfish_neuro)
DimPlot(merfish_neuro, reduction = "umap", group.by = "class_label", pt.size = 0.1, raster = FALSE)

DimPlot(merfish_neuro, reduction = "umap", group.by = "subclass_label", cols = cols_subclass, pt.size = 0.1, raster = FALSE)

DimPlot(merfish_neuro, reduction = "umap", group.by = "cell.type_manual", cols = cols_type, pt.size = 0.1, raster = FALSE)

DimPlot(merfish_neuro, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, raster = FALSE)

table(merfish_neuro$predicted.id)
table(merfish_neuro$predicted.id, merfish_neuro$sub.clustering)
DimPlot(merfish_neuro, reduction = "spatial", group.by = "predicted.id", split.by = "predicted.id")
FeaturePlot(merfish_neuro,"prediction.score.Glutamatergic")
FeaturePlot(merfish_neuro, colnames(predictions)[-1])

## Find markers for remaining clusters==============================================================
Idents(merfish_neuro) <- "integrated_snn_res.0.5"
DefaultAssay(merfish_neuro) = "RNA"
markers_0.5 = FindAllMarkers(merfish_neuro,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             
write.table(markers_0.5, file = "outs/merged/markers_0.5_NeuroOnly.txt", quote = F, sep = "\t", col.names = NA)
to.plot = markers_0.5 %>% group_by(cluster) %>% top_n(12,avg_log2FC) %>% as.data.frame()
cols = rev(colorRampPalette(brewer.pal("RdBu",n = 11))(256))
pdf("outs/merged/top10markers_heatmap.pdf", width = 15, height = 10, useDingbats = F)

## Annotate==============================================================
head(merfish_neuro@meta.data)
DimPlot(merfish_neuro, group.by = "integrated_snn_res.0.5")
merfish_neuro$manual.label =  rep(ifelse(merfish_neuro$integrated_snn_res.0.5 %in% c(15), "Excitatory_Thalamus", merfish_neuro$subclass_label))#Reference is missing this cell type ( thalamus neurons) it clusters distinctively, so I'll add it manually as additional cell type
merfish_neuro$manual.label =  rep(ifelse(merfish_neuro$integrated_snn_res.0.5 %in% c(4), "Inhibitory_BrainStem", merfish_neuro$manual.label))#Reference is missing this cell type, it clusters distinctively, so I'll add it manually as additional cell type
merfish_neuro$manual.label =  rep(ifelse(merfish_neuro$integrated_snn_res.0.5 %in% c(24), "Chat", merfish_neuro$manual.label))#Reference is missing this cell type, it clusters distinctively, so I'll add it manually as additional cell type
DimPlot(merfish_neuro, group.by = "manual.label")


# STRIATUM ONLY -  Subset integrated merfish data to striatum neuronal cells only, and re-run--------------------------------

striatal_clusters = c(2,11) #Cluster 11 is doublets of oligos + striatal neurons, but let's include it anyway for now
merfish_striatum = subset(merfish, integrated_snn_res.0.2 %in% striatal_clusters)
DimPlot(merfish_striatum, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

merfish_striatum <-  NormalizeData(merfish_striatum, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_striatum <- ScaleData(merfish_striatum, features = rownames(merfish_striatum), vars.to.regress = "nCount_RNA")

top30 <- head(VariableFeatures(merfish_striatum), 30)
plot1 <- VariableFeaturePlot(merfish_striatum, pt.size = 0.1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = T)
plot2
#ggsave(plot2, file = "outs/var_features.pdf", device = "pdf")

##PCA=====================================================================
merfish_striatum <-  RunPCA(merfish_striatum, npcs = 100, assay = "integrated")
print(merfish_striatum[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish_striatum, ndims = 100)

VizDimLoadings(merfish_striatum, dims = c(1:10), reduction = "pca")

FeatureScatter(merfish_striatum, "PC_1","nCount_RNA")
DimPlot(merfish_striatum, reduction = "pca", group.by = "orig.ident")

DimHeatmap(merfish_striatum, dims = c(1:10), cells = 500, balanced = TRUE)


##UMAP/tSNE================================================================
merfish_striatum <- RunUMAP(merfish_striatum, 
                            dims = 1:10,
                            seed.use = 6,
                            #n.neighbors = 5,
                            #local.connectivity = 15,
                            #assay = "SCT",
                            #n.components = 3
)

cols_subclass = unique(merfish_striatum$subclass_color)
names(cols_subclass) = unique(merfish_striatum$subclass_label)
DimPlot(merfish_striatum, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, raster = FALSE)

DimPlot(merfish_striatum, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

##Clustering (no label transfer here, instead cluster de novo)=======================================================
library(future)
plan("multisession", workers = 14)
plan() # check
options(future.seed=TRUE)
DefaultAssay(merfish_striatum) = "integrated"
merfish_striatum = FindNeighbors(merfish_striatum, assay = "integrated",
                                 dims = 1:10,
                                 reduction = "pca",
                                 #k.param = 7,
                                 #do.plot = T
)
merfish_striatum@meta.data[grep(pattern = "^integrated_snn", x = names(merfish@meta.data))] = NULL
merfish_striatum$seurat_clusters = NULL
merfish_striatum = FindClusters(merfish_striatum, algorithm = 1, resolution = c(0.05, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
DimPlot(merfish_striatum, reduction = "umap", group.by = c("integrated_snn_res.0.3"), label = T, pt.size = 0.1, raster = FALSE)

DimPlot(merfish_striatum, reduction = "umap", group.by = c("orig.ident"), label = T, pt.size = 0.1, raster = FALSE)

#Check markers for D1 and D2 spiny neurons
FeaturePlot(merfish_striatum,c("Drd1", "Adora2a"))
VlnPlot(merfish_striatum,c("Drd1", "Adora2a"), group.by = "integrated_snn_res.0.3")
FeatureScatter(subset(merfish_striatum,integrated_snn_res.0.3 == 7),  "Drd1", "Adora2a")

DimPlot(merfish_striatum, reduction = "spatial", group.by = c("RNA_snn_res.0.4"), label = T, pt.size = 0.1, raster = FALSE)

clustree(merfish_striatum, prefix ="RNA_snn_res.")

## Find markers==============================================================
Idents(merfish_striatum) <- "integrated_snn_res.0.3"
DefaultAssay(merfish_striatum) = "RNA"
markers_0.3 = FindAllMarkers(merfish_striatum,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             
#write.table(markers_0.3, file = "markers_0.3_StriatumOnly.txt", quote = F, sep = "\t", col.names = NA)

FeaturePlot(merfish_striatum, features = read_clip(), ncol = 3, reduction = "umap", pt.size = 0.1)
FeaturePlot(merfish_striatum, reduction = "umap", features = c("Slc17a7"), pt.size = 0.1)

## Annotate==============================================================
head(merfish_striatum@meta.data)
merfish_striatum$class_label = rep("GABAergic", ncol(merfish_striatum))
merfish_striatum$manual.label = rep(ifelse(merfish_striatum$integrated_snn_res.0.3 %in% c(1,2,4,5,7,8,11), "D1 MSN", "D2 MSN"))
merfish_striatum$manual.label = rep(ifelse(merfish_striatum$integrated_snn_res.0.3 %in% c(10), "doublet_Oligo-MSN", merfish_striatum$manual.label))
DimPlot(merfish_striatum, group.by = "manual.label")

#VASCULAR ONLY-------------------------------
vascular_clusters = c(5,7)
merfish_vascular = subset(merfish, integrated_snn_res.0.2 %in% vascular_clusters)
DimPlot(merfish_vascular, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

merfish_vascular <-  NormalizeData(merfish_vascular, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_vascular <- ScaleData(merfish_vascular, features = rownames(merfish_vascular), vars.to.regress = "nCount_RNA", assay = "RNA")

top30 <- head(VariableFeatures(merfish_vascular), 30)
plot1 <- VariableFeaturePlot(merfish_vascular, pt.size = 0.1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = T)
plot2

##PCA=====================================================================
merfish_vascular <-  RunPCA(merfish_vascular, npcs = 100, assay = "integrated")
print(merfish_vascular[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish_vascular, ndims = 100)

VizDimLoadings(merfish_vascular, dims = c(1:10), reduction = "pca")

FeatureScatter(merfish_vascular, "PC_1","nCount_RNA")
DimPlot(merfish_vascular, reduction = "pca", group.by = "orig.ident")

DimHeatmap(merfish_vascular, dims = c(1:10), cells = 500, balanced = TRUE)

##UMAP/tSNE================================================================
merfish_vascular <- RunUMAP(merfish_vascular, 
                            dims = 1:7,
                            seed.use = 6,
                            #n.neighbors = 5,
                            #local.connectivity = 15,
                            #assay = "SCT",
                            #n.components = 3
)
merfish_vascular <- RunTSNE(merfish_vascular, dims = 1:10) 
#colors = c("#b2182b","#fdae61","#2166ac","#92c5de")
cols_subclass = unique(merfish_vascular$subclass_color)
names(cols_subclass) = unique(merfish_vascular$subclass_label)
DimPlot(merfish_vascular, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, raster = FALSE)

DimPlot(merfish_vascular, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

##Clustering (no label transfer here, instead cluster de novo)=======================================================
library(future)
plan("multisession", workers = 14)
plan() # check
options(future.seed=TRUE)
DefaultAssay(merfish_vascular) = "integrated"
merfish_vascular = FindNeighbors(merfish_vascular, assay = "integrated",
                                 dims = 1:7,
                                 reduction = "pca",
                                 #k.param = 7,
                                 #do.plot = T
)
merfish_vascular@meta.data[grep(pattern = "^integrated_snn", x = names(merfish@meta.data))] = NULL
merfish_vascular$seurat_clusters = NULL
merfish_vascular = FindClusters(merfish_vascular, algorithm = 1, resolution = c(0.05, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
DimPlot(merfish_vascular, reduction = "umap", group.by = c("integrated_snn_res.0.2"), label = T, pt.size = 0.1, raster = FALSE)

DimPlot(merfish_vascular, reduction = "umap", group.by = c("orig.ident"), label = T, pt.size = 0.1, raster = FALSE)

FeaturePlot(merfish_vascular,read_clip())
FeaturePlot(merfish_vascular,c("Rgs5", "blank_mean"))
VlnPlot(merfish_vascular,c("Flt1", "Pecam1", "Kcnj8", "Ace2", "Bmp7", "Prdx6", "Cspg5"), group.by = "integrated_snn_res.0.2")

DimPlot(merfish_vascular, reduction = "spatial", group.by = c("RNA_snn_res.0.4"), label = T, pt.size = 0.1, raster = FALSE)

clustree(merfish_vascular, prefix ="RNA_snn_res.")

saveRDS(merfish_vascular, "merfish_vascular.rds")
## Find markers==============================================================
Idents(merfish_vascular) <- "integrated_snn_res.0.2"
DefaultAssay(merfish_vascular) = "RNA"
markers_0.2 = FindAllMarkers(merfish_vascular,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             
#write.table(markers_0.2, file = "outs/merged/markers_0.2_VascularOnly.txt", quote = F, sep = "\t", col.names = NA)

FeaturePlot(merfish_vascular, features = read_clip(), ncol = 3, reduction = "umap", pt.size = 0.1)

## Annotate==============================================================
head(merfish_vascular@meta.data)
merfish_vascular$class_label = rep("Vascular", ncol(merfish_vascular))
merfish_vascular$manual.label = rep(ifelse(merfish_vascular$integrated_snn_res.0.2 %in% c(0,2,5), "Endothelial", "Other"))
merfish_vascular$manual.label = rep(ifelse(merfish_vascular$integrated_snn_res.0.2 %in% c(1), "Pericytes", merfish_vascular$manual.label))
merfish_vascular$manual.label = rep(ifelse(merfish_vascular$integrated_snn_res.0.2 %in% c(3), "VLMC", merfish_vascular$manual.label))
merfish_vascular$manual.label = rep(ifelse(merfish_vascular$integrated_snn_res.0.2 %in% c(4), "VSMC", merfish_vascular$manual.label))
merfish_vascular$manual.label = rep(ifelse(merfish_vascular$integrated_snn_res.0.2 %in% c(6), "doublet_Vascular-Astrocyte", merfish_vascular$manual.label))
DimPlot(merfish_vascular, group.by = "manual.label")

#IMMUNE ONLY-------------------------------
immune_clusters = c(6)
merfish_immune = subset(merfish, integrated_snn_res.0.2 %in% immune_clusters)
DimPlot(merfish_immune, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, raster = FALSE)

merfish_immune <-  NormalizeData(merfish_immune, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_immune <- ScaleData(merfish_immune, features = rownames(merfish_immune), vars.to.regress = "nCount_RNA", assay = "RNA")

##Since Immune cells are most contaminated by RNA from other cells, let's process them again---------------
###Prepare for regressing out gene marker scores of other cell types===========================================
DEgenes_overlap = list()
DEgenes = list()
markers_0.2 = read.csv("outs/merged/markers_0.2.txt", sep = "\t", row.names = 1)
DimPlot(merfish)
#subset merfish to only relevant major clusters
merfish_selected = subset(merfish, integrated_snn_res.0.2 %in% c(0,1,2,3,4,5,6,7,8,9,12,14))
Idents(merfish_selected) = "integrated_snn_res.0.2"
DimPlot(merfish_selected)
for( i in setdiff(as.numeric(levels(Idents(merfish_selected))), immune_clusters)){
  temp = markers_0.2 %>% filter(cluster ==i, avg_log2FC >1 & p_val_adj < 0.01) #use markers from full set of all cells
  name = paste0("cluster_",i)
  DEgenes[[name]]= FindMarkers(merfish_selected, ident.1 = i, ident.2 = immune_clusters,
                               logfc.threshold = 1,
                               test.use = "wilcox",
                               min.diff.pct = 0.3,
                               only.pos = T)
  DEgenes_overlap[[name]] = intersect(rownames(DEgenes[[name]]), temp$gene)

}

merfish_immune = AddModuleScore(merfish_immune, DEgenes_overlap, name = "set_", seed = 42, search = F, assay = "RNA", ctrl = 5)
set_names <- names(merfish_immune@meta.data)[grep(pattern = "^set_", x = names(merfish_immune@meta.data))]
new_meta = merfish_immune@meta.data[set_names]
colnames(new_meta) = names(DEgenes_overlap)
merfish_immune@meta.data[set_names] = NULL
merfish_immune = AddMetaData(merfish_immune, metadata = new_meta)

#Plot cell marker scores
FeaturePlot(merfish_immune, names(DEgenes_overlap), raster = TRUE)

#Re-integrate immune cells over 3 sections
##Integration---------------------------------------------------
# split the dataset into a list of seurat objects (each treatment will be separate object)
options(future.globals.maxSize= 2891289600)
merfish_immune.list <- SplitObject(merfish_immune, split.by = "orig.ident")
variables.to.regress = names(DEgenes_overlap)[which(names(DEgenes_overlap) != "cluster_7")] #Exclude cluster_7 as it is essentially same as cluster_5 (vascular cells)
# normalize and scale and regress variables for each library separately
for (i in 1:length(merfish_immune.list)) {
  merfish_immune.list[[i]] <- NormalizeData(merfish_immune.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  merfish_immune.list[[i]] <- ScaleData(merfish_immune.list[[i]], features = rownames(merfish_immune.list[[i]]), vars.to.regress = c(variables.to.regress, "nCount_RNA"))
}

gene.annot = read.table("data/gene.annotation.txt", header = T)
rownames(gene.annot) = gene.annot$Gene
Microglial = gene.annot[gene.annot$annotation == "Microglial", "Gene"]
nonMicroglial = gene.annot[gene.annot$annotation == "nonMicroglial", "Gene"]
merfish_immune.features = Microglial

#Find integration anchors
merfish_immune.anchors <- FindIntegrationAnchors(object.list = merfish_immune.list,
                                                 normalization.method = "LogNormalize",
                                                 assay = c("RNA", "RNA", "RNA"),
                                                 anchor.features = merfish_immune.features,
                                                 reduction = "cca",
                                                 dims = 1:30, k.anchor = 5) 
# Integrate data, command creates an 'integrated' data assay
all_features <- lapply(merfish_immune.list, row.names)
all_features <- Reduce(intersect, all_features)
merfish_immune <- IntegrateData(anchorset = merfish_immune.anchors,
                                normalization.method = "LogNormalize",
                                features.to.integrate = all_features,
                                dims = 1:30
)

merfish_immune <-  NormalizeData(merfish_immune, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_immune <- ScaleData(merfish_immune, features = rownames(merfish_immune), assay = "RNA", vars.to.regress = c("nCount_RNA", "orig.ident"))


##PCA=====================================================================
DefaultAssay(merfish_immune) = "integrated"
merfish_immune = ScaleData(merfish_immune, features = rownames(merfish_immune), vars.to.regress = c(variables.to.regress, "nCount_RNA"))

selected.genes = c("Selplg", "Laptm5", "Csf1r", "Tmem119", "C1qa", "C1qb",
                   "Mrc1", "Cd163",
                   "Cd3e", "Cd28", "Cd247", "Tcf7", "Cd4", "Cd8a",
                   "Stat1", "Rsad2", "Ifng", "Usp18", "Ifit1")

merfish_immune <-  RunPCA(merfish_immune, npcs = 35, features =merfish_immune.features,  assay = "integrated")
print(merfish_immune[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish_immune, ndims = 100)

VizDimLoadings(merfish_immune, dims = c(1:10), reduction = "pca")

FeatureScatter(merfish_immune, "PC_1","nCount_RNA")
DimPlot(merfish_immune, reduction = "pca", group.by = "orig.ident")

DimHeatmap(merfish_immune, dims = c(1:10), cells = 500, balanced = TRUE)

##UMAP/tSNE================================================================
merfish_immune <- RunUMAP(merfish_immune,
                          dims = 1:20,
                          seed.use = 6,
                          #n.neighbors = 5,
                          #local.connectivity = 15,
                          #assay = "SCT",
                          #n.components = 3
)
merfish_immune <- RunTSNE(merfish_immune, dims = 1:10) 
cols_subclass = unique(merfish_immune$subclass_color)
names(cols_subclass) = unique(merfish_immune$subclass_label)
DimPlot(merfish_immune, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE)

DimPlot(merfish_immune, reduction = "umap", group.by = "hemisphere", pt.size = 0.1, raster = FALSE)

#Check expression of some marker genes
DefaultAssay(merfish_immune) = "RNA"
FeaturePlot(merfish,c("Cd3e", "Cd4", "Cd8a", "Tcf7", "Cd247", "Cd28", "Ifng"), order = T)
FeaturePlot(merfish_immune,c("Clec7a", "Itgax", "Apoe", "Cd74", "Spp1"))
FeaturePlot(merfish_immune,c("Gpnmb", "Lgals3", "Lpl"))
FeaturePlot(merfish_immune,c("Plin2", "Plin3", "Soat1", "Anxa5"))
FeaturePlot(merfish_immune,c("Mrc1", "Cd163", "Msr1"))
FeaturePlot(merfish_immune,c("Stat1", "Rsad2", "Usp18", "Ifit1"))
FeaturePlot(merfish_immune,c("Slc17a7", "Plp1", "Gad2", "Flt1", "Drd1"))
FeaturePlot(merfish_immune,variables.to.regress)

VlnPlot(merfish,c("Cd3e", "Cd4", "Cd8a", "Tcf7", "Cd247", "Cd28", "Ifng"), pt.size = 0)

##Clustering (no label transfer here, instead cluster de novo)=======================================================
library(future)
plan("multisession", workers = 14)
plan() # check
options(future.seed=TRUE)
DefaultAssay(merfish_immune) = "integrated"
merfish_immune = FindNeighbors(merfish_immune, assay = "integrated",
                               dims = 1:20,
                               reduction = "pca",
                               k.param = 7,
                               #do.plot = T
)
merfish_immune@meta.data[grep(pattern = "^integrated_snn", x = names(merfish@meta.data))] = NULL
merfish_immune$seurat_clusters = NULL
merfish_immune = FindClusters(merfish_immune, algorithm = 1, resolution = c(seq(0.1, 1.6, by = 0.1)))
DimPlot(merfish_immune, reduction = "umap", group.by = c("integrated_snn_res.0.2"), label = T, pt.size = 0.5, raster = FALSE)

Tcell_markers = c("Cd3e", "Cd4", "Cd8a", "Tcf7", "Cd247", "Cd28", "Ifng")
Macro_markers = c("Mrc1", "Cd163")
scores = list(Tcell = Tcell_markers, Macro = Macro_markers)

DefaultAssay(merfish_immune) = "RNA"
merfish_immune = AddModuleScore(merfish_immune, scores, name = "set_", seed = 42, search = F, assay = "RNA", ctrl = 5)
set_names <- names(merfish_immune@meta.data)[grep(pattern = "^set_", x = names(merfish_immune@meta.data))]
new_meta = merfish_immune@meta.data[set_names]
colnames(new_meta) = names(scores)
merfish_immune@meta.data[set_names] = NULL
merfish_immune = AddMetaData(merfish_immune, metadata = new_meta)
#merfish_immune@meta.data[names(scores)] = NULL
FeaturePlot(merfish_immune,names(scores))
FeatureScatter(merfish_immune, "Selplg", "Mrc1")

clustree(merfish_immune, prefix ="integrated_snn_res.")

saveRDS(merfish_immune, "data/merfish_immune.rds")
## Find markers==============================================================
Idents(merfish_immune) <- "integrated_snn_res.1.5"
DefaultAssay(merfish_immune) = "RNA"
markers_1.5 = FindAllMarkers(merfish_immune,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             

## Annotate=============================================================
head(merfish_immune@meta.data)

merfish_immune$class_label = rep("Immune", ncol(merfish_immune))
merfish_immune$manual.label = rep(ifelse(merfish_immune$integrated_snn_res.1.5 %in% c(20), "Macrophages", "Microglia"))
merfish_immune$manual.label = rep(ifelse(merfish_immune$integrated_snn_res.1.5 %in% c(21), "T-cells", merfish_immune$manual.label))
DimPlot(merfish_immune, group.by = "manual.label", pt.size = 0.5)
VlnPlot(merfish_immune, selected.genes, group.by = "manual.label", pt.size = 0.1)

#MICROGLIA ONLY-------------------------------
Idents(merfish_immune) = "integrated_snn_res.1.5"
non_microglia_clusters = c(20, 21)
merfish_microglia = subset(merfish_immune, integrated_snn_res.1.5 %in% non_microglia_clusters, invert = TRUE)
DimPlot(merfish_microglia, reduction = "umap", group.by = "integrated_snn_res.1.5", pt.size = 0.1, raster = FALSE)

merfish_microglia <-  NormalizeData(merfish_microglia, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_microglia <- ScaleData(merfish_microglia, features = rownames(merfish_microglia), vars.to.regress = "nCount_RNA", assay = "RNA")

##reprocess microglia ---------------
###Prepare for regressing out gene marker scores of other cell types===========================================
merfish_microglia = AddModuleScore(merfish_microglia, DEgenes_overlap, name = "set_", seed = 42, search = F, assay = "RNA", ctrl = 5)
set_names <- names(merfish_microglia@meta.data)[grep(pattern = "^set_", x = names(merfish_microglia@meta.data))]
new_meta = merfish_microglia@meta.data[set_names]
colnames(new_meta) = names(DEgenes_overlap)
merfish_microglia@meta.data[set_names] = NULL
merfish_microglia = AddMetaData(merfish_microglia, metadata = new_meta)

FeaturePlot(merfish_microglia, names(DEgenes_overlap), raster = TRUE) 

#Re-integrate microglia cells over 3 sections
##Integration---------------------------------------------------
# split the dataset into a list of seurat objects (each treatment will be separate object)
options(future.globals.maxSize= 2891289600)
merfish_microglia.list <- SplitObject(merfish_microglia, split.by = "orig.ident")
variables.to.regress = names(DEgenes_overlap)[which(names(DEgenes_overlap) != "cluster_7")] #Exclude cluster_7 as it is essentially same as cluster_5
# normalize and scale and regress variables for each library separately
for (i in 1:length(merfish_microglia.list)) {
  merfish_microglia.list[[i]] <- NormalizeData(merfish_microglia.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  merfish_microglia.list[[i]] <- ScaleData(merfish_microglia.list[[i]], features = rownames(merfish_microglia.list[[i]]), vars.to.regress = c(variables.to.regress, "nCount_RNA"))
}


Microglial = gene.annot[gene.annot$annotation == "Microglial", "Gene"]
nonMicroglial = gene.annot[gene.annot$annotation == "nonMicroglial", "Gene"]
merfish_microglia.features = Microglial

for (name in names(merfish_microglia.list)) { 
  print(nrow(merfish_microglia.list[[name]]@assays[["RNA"]]@scale.data))#check the number of genes in "scale.data" slot of "RNA" assay
  print(length(VariableFeatures(merfish_microglia.list[[name]])))#check the number of variable features for each dataset
}

merfish_microglia.anchors <- FindIntegrationAnchors(object.list = merfish_microglia.list,
                                                    normalization.method = "LogNormalize",
                                                    assay = c("RNA", "RNA", "RNA"),
                                                    anchor.features = merfish_microglia.features,
                                                    reduction = "cca",
                                                    dims = 1:30, k.anchor = 5) 
# Integrate data, command creates an 'integrated' data assay
all_features <- lapply(merfish_microglia.list, row.names)
all_features <- Reduce(intersect, all_features)
merfish_microglia <- IntegrateData(anchorset = merfish_microglia.anchors,
                                   normalization.method = "LogNormalize",
                                   features.to.integrate = all_features,
                                   dims = 1:30
)

merfish_microglia <-  NormalizeData(merfish_microglia, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish_microglia <- ScaleData(merfish_microglia, features = rownames(merfish_microglia), assay = "RNA", vars.to.regress = c("nCount_RNA", "orig.ident"))

##PCA=====================================================================
DefaultAssay(merfish_microglia) = "integrated"
merfish_microglia = ScaleData(merfish_microglia, features = rownames(merfish_microglia), vars.to.regress = c(variables.to.regress, "nCount_RNA"))

selected.genes = c("Selplg", "Laptm5", "Csf1r", "Tmem119", "C1qa", "C1qb",
                   "Mrc1", "Cd163",
                   "Cd3e", "Cd28", "Cd247", "Tcf7", "Cd4", "Cd8a",
                   "Stat1", "Rsad2", "Ifng", "Usp18", "Ifit1")

merfish_microglia <-  RunPCA(merfish_microglia, npcs = 35, features =merfish_microglia.features,  assay = "integrated")
print(merfish_microglia[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(merfish_microglia, ndims = 100)

VizDimLoadings(merfish_microglia, dims = c(1:10), reduction = "pca")

FeatureScatter(merfish_microglia, "PC_1","nCount_RNA")
DimPlot(merfish_microglia, reduction = "pca", group.by = "orig.ident")

DimHeatmap(merfish_microglia, dims = c(1:10), cells = 500, balanced = TRUE)

##UMAP/tSNE================================================================
merfish_microglia <- RunUMAP(merfish_microglia,
                             dims = 1:20,
                             seed.use = 6,
                             #n.neighbors = 5,
                             #local.connectivity = 15,
                             #assay = "SCT",
                             #n.components = 3
)
DimPlot(merfish_microglia, reduction = "umap", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)

DimPlot(merfish_microglia, reduction = "umap", group.by = "region.coarse", pt.size = 0.5, raster = FALSE)

DefaultAssay(merfish_microglia) = "RNA"
FeaturePlot(merfish_microglia,c("Cd3e", "Cd4", "Cd8a", "Tcf7", "Cd247", "Cd28", "Ifng"), order = T)
FeaturePlot(merfish_microglia,c("Clec7a", "Itgax", "Apoe", "Cd74", "Spp1"))
FeaturePlot(merfish_microglia,c("Gpnmb", "Lgals3", "Lpl"))
FeaturePlot(merfish_microglia,c("Plin2", "Plin3", "Soat1", "Abca1", "Abcg1"))
FeaturePlot(merfish_microglia,c("Mrc1", "Cd163", "Msr1"))
FeaturePlot(merfish_microglia,c("Stat1", "Rsad2", "Usp18", "Ifit1"))
FeaturePlot(merfish_microglia,c("Slc17a7", "Plp1", "Gad2", "Flt1", "Drd1", "nCount_RNA"))
FeaturePlot(merfish_microglia,variables.to.regress)

##Clustering (no label transfer here, instead cluster de novo)=======================================================
library(future)
plan("multisession", workers = 14)
plan() # check
options(future.seed=TRUE)
DefaultAssay(merfish_microglia) = "integrated"
merfish_microglia = FindNeighbors(merfish_microglia, assay = "integrated",
                                  dims = 1:20,
                                  reduction = "pca",
                                  #k.param = 7,
                                  #do.plot = T
)
merfish_microglia@meta.data[grep(pattern = "^integrated_snn", x = names(merfish@meta.data))] = NULL
merfish_microglia$seurat_clusters = NULL
merfish_microglia = FindClusters(merfish_microglia, algorithm = 1, resolution = c(seq(0.1, 1.6, by = 0.1)))
DimPlot(merfish_microglia, reduction = "umap", group.by = c("integrated_snn_res.0.3"), label = T, pt.size = 0.5, raster = FALSE)


saveRDS(merfish_microglia, "data/merfish_microglia.rds")

## Find markers==============================================================
Idents(merfish_microglia) <- "integrated_snn_res.0.3"
DefaultAssay(merfish_microglia) = "RNA"
markers_0.3 = FindAllMarkers(merfish_microglia,
                             logfc.threshold = 0.3,
                             test.use = "wilcox",
                             min.pct = 0.2,
                             min.diff.pct = -Inf,
                             only.pos = T,
                             future.seed = TRUE
)                             
write.table(markers_0.3, file = "markers_0.3_MicrogliaOnly-reintegrated.txt", quote = F, sep = "\t", col.names = NA)

## Annotate=============================================================
head(merfish_microglia@meta.data)

DimPlot(merfish_microglia, group.by = "integrated_snn_res.1.5", pt.size = 0.5, cells.highlight = colnames(subset(merfish_microglia,integrated_snn_res.1.5 == 16 )))
FeaturePlot(merfish_microglia, "Gpnmb", order = T, pt.size = 0.6)

merfish_microglia$manual.label = rep(ifelse(merfish_microglia$integrated_snn_res.0.3 %in% c(0,3,4,6), "Micro_Homeostatic", "Micro_Activated_1"))
merfish_microglia$manual.label = rep(ifelse(merfish_microglia$integrated_snn_res.0.3 %in% c(5), "Micro_Interferon", merfish_microglia$manual.label))
merfish_microglia$manual.label = rep(ifelse(merfish_microglia$integrated_snn_res.1.5 %in% c(16), "Micro_Activated_2", merfish_microglia$manual.label))

DimPlot(merfish_microglia, group.by = "manual.label", pt.size = 0.5)

#FINISH ANNOTATION----------------------------------------------------
##Finish neuronal labeling--------------
table(ABA_metadata$class_label, ABA_metadata$subclass_label)
GABAergic = ABA_metadata %>% filter(class_label == "GABAergic") %>% pull(subclass_label) %>% unique()
GABAergic = c(GABAergic, "Chat", "Inhibitory_BrainStem")
Glutamatergic = ABA_metadata %>% filter(class_label == "Glutamatergic") %>% pull(subclass_label) %>% unique()
Glutamatergic = c(Glutamatergic, "Excitatory_Thalamus")

merfish_neuro$cell.class_manual = rep(ifelse(merfish_neuro$manual.label %in% GABAergic , "GABAergic", "Glutamatergic"))
DimPlot(merfish_neuro, group.by = "cell.class_manual", pt.size = 0.5)

merfish_neuro$cell.type_manual = merfish_neuro$manual.label
DimPlot(merfish_neuro, group.by = "cell.type_manual", pt.size = 0.5)

merfish_neuro$cell.cluster_manual = merfish_neuro$cell.type_manual
DimPlot(merfish_neuro, group.by = "cell.cluster_manual", pt.size = 0.5)


##Finish striatum labeling--------------
table(merfish_striatum$manual.label)
merfish_striatum$cell.class_manual = rep(ifelse(merfish_striatum$manual.label %in% c("D1 MSN", "D2 MSN"), "GABAergic", merfish_striatum$manual.label))
DimPlot(merfish_striatum, group.by = "cell.class_manual", pt.size = 0.5)

merfish_striatum$cell.type_manual = merfish_striatum$manual.label
DimPlot(merfish_striatum, group.by = "cell.type_manual", pt.size = 0.5)

merfish_striatum$cell.cluster_manual = merfish_striatum$cell.type_manual
DimPlot(merfish_striatum, group.by = "cell.cluster_manual", pt.size = 0.5)

##Finish vascular labeling--------------
table(merfish_vascular$manual.label)
merfish_vascular$cell.class_manual = rep(ifelse(merfish_vascular$manual.label != "doublet_Vascular-Astrocyte", "Vascular", merfish_vascular$manual.label))
DimPlot(merfish_vascular, group.by = "cell.class_manual", pt.size = 0.5)

merfish_vascular$cell.type_manual = merfish_vascular$manual.label
DimPlot(merfish_vascular, group.by = "cell.type_manual", pt.size = 0.5)

merfish_vascular$cell.cluster_manual = merfish_vascular$cell.type_manual
DimPlot(merfish_vascular, group.by = "cell.cluster_manual", pt.size = 0.5)

##Finish immune labeling--------------
table(merfish_immune$manual.label)
merfish_immune$cell.class_manual = rep("Immune", nrow(merfish_immune@meta.data))
DimPlot(merfish_immune, group.by = "cell.class_manual", pt.size = 0.5)

merfish_immune$cell.type_manual = merfish_immune$manual.label
DimPlot(merfish_immune, group.by = "cell.type_manual", pt.size = 0.5)

merfish_microglia$cell.class_manual = rep("Immune", nrow(merfish_microglia@meta.data))
DimPlot(merfish_microglia, group.by = "cell.class_manual", pt.size = 0.5)

merfish_microglia$cell.type_manual = rep("Microglia", nrow(merfish_microglia@meta.data))
DimPlot(merfish_microglia, group.by = "cell.type_manual", pt.size = 0.5)

merfish_microglia$cell.cluster_manual = merfish_microglia$manual.label
DimPlot(merfish_microglia, group.by = "cell.cluster_manual", pt.size = 0.5)

temp = merfish_microglia@meta.data %>% select("manual.label") %>% rename(temp = "manual.label")
head(temp)
merfish_immune = AddMetaData(merfish_immune, metadata = temp)
DimPlot(merfish_immune, group.by = "temp", pt.size = 0.5)

merfish_immune$cell.cluster_manual = rep(ifelse(is.na(merfish_immune$temp), merfish_immune$cell.type_manual, merfish_immune$temp))
DimPlot(merfish_immune, group.by = "cell.cluster_manual", pt.size = 0.5)
merfish_immune$temp = NULL
table(merfish_immune$cell.cluster_manual, merfish_immune$cell.type_manual)

##Finish remaining cell type labeling--------------
#REMOVE Cluster 10 from all cells (these are low quality cells originating from crack in sample b3s20)
merfish = subset(merfish,integrated_snn_res.0.2 == 10, invert = TRUE)
DimPlot(merfish, group.by = "cell.type_manual", pt.size = 0.1, raster = FALSE, label = TRUE)

#Collect annotation meta from sub objects
cell.annotation = merfish_neuro@meta.data %>% select(c("cell.class_manual", "cell.type_manual", "cell.cluster_manual"))
cell.annotation = cell.annotation %>% rbind(.,merfish_striatum@meta.data %>% select(c("cell.class_manual", "cell.type_manual", "cell.cluster_manual")))
cell.annotation = cell.annotation %>% rbind(.,merfish_vascular@meta.data %>% select(c("cell.class_manual", "cell.type_manual", "cell.cluster_manual")))
cell.annotation = cell.annotation %>% rbind(.,merfish_immune@meta.data %>% select(c("cell.class_manual", "cell.type_manual", "cell.cluster_manual")))
head(cell.annotation)

merfish = AddMetaData(merfish, metadata = cell.annotation)
cell.annotation2 = merfish@meta.data %>% filter(sample.name %notin% rownames(cell.annotation)) %>% select(integrated_snn_res.0.2)
cell.annotation2$cell.class_manual = rep(ifelse(cell.annotation2$integrated_snn_res.0.2 == 3, "Astrocytes", "notDefined"))
cell.annotation2$cell.class_manual = rep(ifelse(cell.annotation2$integrated_snn_res.0.2 == 9, "OPCs", cell.annotation2$cell.class_manual))
cell.annotation2$cell.class_manual = rep(ifelse(cell.annotation2$integrated_snn_res.0.2 == 12, "Ependymal", cell.annotation2$cell.class_manual))
cell.annotation2$cell.class_manual = rep(ifelse(cell.annotation2$integrated_snn_res.0.2 == 15, "NeuralPrecursors", cell.annotation2$cell.class_manual))

cell.annotation2$cell.type_manual = cell.annotation2$cell.class_manual
cell.annotation2$cell.cluster_manual = cell.annotation2$cell.class_manual
cell.annotation2$integrated_snn_res.0.2 = NULL
head(cell.annotation2)

cell.annotation = rbind(cell.annotation, cell.annotation2)
merfish = AddMetaData(merfish, metadata = cell.annotation)

DimPlot(merfish, group.by = "cell.cluster_manual", pt.size = 0.1, label = TRUE)
