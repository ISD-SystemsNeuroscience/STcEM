library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratData)
library(plotly)
library(crosstalk)
library(DT)
library(data.table)

#Read in data-------------------------------------
setwd(".")
b3s20_matrix = read.csv("data/b3s20/cell_by_gene.csv")
b3s20_matrix[1:10, 1:10]
colnames(b3s20_matrix)
b3s20_matrix$cellID = paste0("b3s20_", b3s20_matrix$X)
length(unique(b3s20_matrix$cellID)) == nrow(b3s20_matrix) #check if all cells have unique name
rownames(b3s20_matrix) = b3s20_matrix$cellID
b3s20_matrix$X = NULL
b3s20_matrix$cellID = NULL

b3s20_meta = read.csv("data/b3s20/cell_metadata.csv")
b3s20_meta[1:20, 1:ncol(b3s20_meta)]
b3s20_meta$cellID = paste0("b3s20_", b3s20_meta$X)
length(unique(b3s20_meta$cellID)) == nrow(b3s20_meta) #check if all cells have unique name
rownames(b3s20_meta) = b3s20_meta$cellID
b3s20_meta$X = NULL

table(rownames(b3s20_matrix) == rownames(b3s20_meta))
table(rownames(b3s20_matrix) %in% rownames(b3s20_meta))
b3s20_matrix = b3s20_matrix[order(match(rownames(b3s20_matrix), rownames(b3s20_meta))),] #reorder rows of matrix to match metadata

b3s20_transcripts = fread("data/b3s20/detected_transcripts.csv", data.table = FALSE) #large file, use data.table::fread() to read it fast
b3s20_transcripts[1:100, 1:10]
b3s20_transcripts$uniqueID = paste0("b2s20_", rep(1:nrow(b3s20_transcripts)))
rownames(b3s20_transcripts) = b3s20_transcripts$uniqueID
b3s20_transcripts$V1 = NULL

#Subsample transcripts for visualization purposes
theme_set(theme_bw())
transcripts_sampled = sample(rownames(b3s20_transcripts), 100000)
transcripts_sampled = b3s20_transcripts[transcripts_sampled,]
p1 = ggplot()+geom_point(data=transcripts_sampled, aes(x = global_x, y = global_y, color = gene), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p1
#ggsave(plot = p1, filename = "outs/b3s20/transcripts_downsampled.pdf", height = 7, width = 9)

#Plotting transcript and cell positions
p2 = ggplot()+geom_point(data=b3s20_meta, aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p2
#ggsave(plot = p2, filename = "outs/b3s20/cells.pdf", height = 7, width = 9)

# Transform (rotate) coordinates of cells, transcripts and images-------------------------
## Rotate cells
embeddings = as.matrix(b3s20_meta[,c("center_x", "center_y")])
ggplot(as.data.frame(embeddings)) + geom_point(aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)
embeddings.rotated <- as.data.frame(spdep::Rotation(embeddings, 228.5*pi/180))
colnames(embeddings.rotated) = c("rotated_x", "rotated_y")
ggplot(embeddings.rotated) + geom_point(aes(x = rotated_x, y = rotated_y), size = 0.1, alpha = 0.4)
#ggsave(filename = "outs/b3s20/cells_rotated.pdf", height = 7, width = 9)

head(embeddings.rotated)
head(b3s20_meta)
embeddings.rotated$cellID = rownames(embeddings.rotated)
b3s20_meta = left_join(b3s20_meta, embeddings.rotated, by = "cellID")
rownames(b3s20_meta) = b3s20_meta$cellID

## Rotate ALL transcripts !--------------------------------------------------
transcripts.rotated = as.matrix(b3s20_transcripts[,c("global_x", "global_y")])
transcripts.rotated = as.data.frame(spdep::Rotation(transcripts.rotated, 228.5*pi/180))
colnames(transcripts.rotated) = c("rotated_x", "rotated_y")
#Add rotated positions to original transcript files
transcripts.rotated$uniqueID = rownames(transcripts.rotated)
b3s20_transcripts = left_join(b3s20_transcripts, transcripts.rotated, by = "uniqueID")
head(b3s20_transcripts)
rownames(b3s20_transcripts) = b3s20_transcripts$uniqueID

#QC and exploratory-------------------------------------
##calculate blank signals
blank_probes = grepl(pattern = "Blank*",colnames(b3s20_matrix))
table(blank_probes)
table(rownames(b3s20_meta) == rownames(b3s20_matrix)) #Must be all TRUE
b3s20_meta$blank_sum = rowSums(b3s20_matrix[blank_probes])
b3s20_meta$blank_mean = rowMeans(b3s20_matrix[blank_probes])

##create Seurat
merfish = CreateSeuratObject(counts = t(b3s20_matrix[!blank_probes]), project = "b3s20", meta.data = b3s20_meta)
p7 = VlnPlot(merfish, features = c("nCount_RNA","nFeature_RNA", "volume", "blank_sum", "blank_mean"), pt.size = 0)
p7
#ggsave(plot = p7, filename = "outs/b3s20/QC_violin.pdf", height = 7, width = 9)

# Add physical spatial dimension to seurat object as dimensional reduction embedding
spatial.dim = as.matrix(merfish@meta.data[c("rotated_x", "rotated_y")])
colnames(spatial.dim) = c("spatial_1", "spatial_2")
spatial = CreateDimReducObject(embeddings = spatial.dim, key = "spatial_")
merfish@reductions[["spatial"]] = spatial

DimPlot(merfish, reduction = "spatial", raster = FALSE)
DimPlot(merfish, reduction = "spatial",raster = FALSE)

#QC plots
cols = viridisLite::turbo(n = 256)
p8 = FeaturePlot(merfish, c("nCount_RNA", "nFeature_RNA","volume", "blank_sum", "blank_mean"),
            reduction = "spatial", pt.size = 0.1, min.cutoff = "q1", max.cutoff = "q99")& scale_colour_gradientn(colors = cols)
p8

FeatureScatter(merfish, "nCount_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "nCount_RNA")
FeatureScatter(merfish, "nFeature_RNA", "blank_mean")

ggplot(merfish@meta.data, aes(log(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(log2(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log2(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(nCount_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nCount_RNA, c(0.01,0.99)))
quantile(merfish$nCount_RNA, c(0.01,0.99))
median(merfish$nCount_RNA)

ggplot(merfish@meta.data, aes(nFeature_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nFeature_RNA, c(0.01,0.99)))
quantile(merfish$nFeature_RNA, c(0.01,0.99))
median(merfish$nFeature_RNA)

ggplot(merfish@meta.data, aes(blank_mean))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$blank_mean, c(0.01,0.99)))
quantile(merfish$blank_mean, c(0.01,0.99))
median(merfish$blank_mean)

#create QC column labeling good and bad cells
good.cells = merfish@meta.data %>% filter(volume <2500 & volume >=40 & nCount_RNA <= 2500 & nCount_RNA >= 30 & nFeature_RNA >= 5 & blank_mean <= 1)
merfish$QCpass = ifelse(colnames(merfish) %in% rownames(good.cells), "TRUE", "FALSE")
table(merfish$QCpass)


## Data subsetting, normalization, scaling------------
merfish = subset(merfish, cells = rownames(good.cells))

merfish <-  NormalizeData(merfish, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish <- ScaleData(merfish, features = rownames(merfish), assay = "RNA", vars.to.regress = "nCount_RNA")

#Save Seurat object
saveRDS(merfish, "data/merfish_b3s20.rds")

#======================Repeat for section b2s20================================

#Read in data-------------------------------------
setwd("~/GOKCE LAB/Code/MERFISH_LPC/github/")
b2s20_matrix = read.csv("data/b2s20/cell_by_gene.csv")
b2s20_matrix[1:10, 1:10]
colnames(b2s20_matrix)
b2s20_matrix$cellID = paste0("b2s20_", b2s20_matrix$X)
length(unique(b2s20_matrix$cellID)) == nrow(b2s20_matrix) #check if all cells have unique name
rownames(b2s20_matrix) = b2s20_matrix$cellID
b2s20_matrix$X = NULL
b2s20_matrix$cellID = NULL

b2s20_meta = read.csv("data/b2s20/cell_metadata.csv")
b2s20_meta[1:20, 1:ncol(b2s20_meta)]
b2s20_meta$cellID = paste0("b2s20_", b2s20_meta$X)
length(unique(b2s20_meta$cellID)) == nrow(b2s20_meta) #check if all cells have unique name
rownames(b2s20_meta) = b2s20_meta$cellID
b2s20_meta$X = NULL

table(rownames(b2s20_matrix) == rownames(b2s20_meta))
table(rownames(b2s20_matrix) %in% rownames(b2s20_meta))
b2s20_matrix = b2s20_matrix[order(match(rownames(b2s20_matrix), rownames(b2s20_meta))),] #reorder rows of matrix to match metadata

b2s20_transcripts = fread("data/b2s20/detected_transcripts.csv", data.table = FALSE) #large file, use data.table::fread() to read it fast
b2s20_transcripts[1:100, 1:10]
b2s20_transcripts$uniqueID = paste0("b2s20_", rep(1:nrow(b2s20_transcripts)))
rownames(b2s20_transcripts) = b2s20_transcripts$uniqueID
b2s20_transcripts$V1 = NULL

#Subsample transcripts for visualization purposes
theme_set(theme_bw())
transcripts_sampled = sample(rownames(b2s20_transcripts), 100000)
transcripts_sampled = b2s20_transcripts[transcripts_sampled,]
p1 = ggplot()+geom_point(data=transcripts_sampled, aes(x = global_x, y = global_y, color = gene), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p1
#ggsave(plot = p1, filename = "outs/b2s20/transcripts_downsampled.pdf", height = 7, width = 9)

#Plotting transcript and cell positions
p2 = ggplot()+geom_point(data=b2s20_meta, aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p2
#ggsave(plot = p2, filename = "outs/b2s20/cells.pdf", height = 7, width = 9)

# Transform (rotate) coordinates of cells, transcripts and images-------------------------
## Rotate cells
embeddings = as.matrix(b2s20_meta[,c("center_x", "center_y")])
ggplot(as.data.frame(embeddings)) + geom_point(aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)
embeddings.rotated <- as.data.frame(spdep::Rotation(embeddings, 170*pi/180))
colnames(embeddings.rotated) = c("rotated_x", "rotated_y")
ggplot(embeddings.rotated) + geom_point(aes(x = rotated_x, y = rotated_y), size = 0.1, alpha = 0.4)
#ggsave(filename = "outs/b2s20/cells_rotated.pdf", height = 7, width = 9)

head(embeddings.rotated)
head(b2s20_meta)
embeddings.rotated$cellID = rownames(embeddings.rotated)
b2s20_meta = left_join(b2s20_meta, embeddings.rotated, by = "cellID")
rownames(b2s20_meta) = b2s20_meta$cellID

## Rotate ALL transcripts !--------------------------------------------------
transcripts.rotated = as.matrix(b2s20_transcripts[,c("global_x", "global_y")])
transcripts.rotated = as.data.frame(spdep::Rotation(transcripts.rotated, 170*pi/180))
colnames(transcripts.rotated) = c("rotated_x", "rotated_y")
#Add rotated positions to original transcript files
transcripts.rotated$uniqueID = rownames(transcripts.rotated)
b2s20_transcripts = left_join(b2s20_transcripts, transcripts.rotated, by = "uniqueID")
head(b2s20_transcripts)
rownames(b2s20_transcripts) = b2s20_transcripts$uniqueID

#QC and exploratory-------------------------------------
##calculate blank signals
blank_probes = grepl(pattern = "Blank*",colnames(b2s20_matrix))
table(blank_probes)
table(rownames(b2s20_meta) == rownames(b2s20_matrix)) #Must be all TRUE
b2s20_meta$blank_sum = rowSums(b2s20_matrix[blank_probes])
b2s20_meta$blank_mean = rowMeans(b2s20_matrix[blank_probes])

##create Seurat
merfish = CreateSeuratObject(counts = t(b2s20_matrix[!blank_probes]), project = "b2s20", meta.data = b2s20_meta)
p7 = VlnPlot(merfish, features = c("nCount_RNA","nFeature_RNA", "volume", "blank_sum", "blank_mean"), pt.size = 0)
p7
#ggsave(plot = p7, filename = "outs/b2s20/QC_violin.pdf", height = 7, width = 9)

# Add physical spatial dimension to seurat object as dimensional reduction embedding
spatial.dim = as.matrix(merfish@meta.data[c("rotated_x", "rotated_y")])
colnames(spatial.dim) = c("spatial_1", "spatial_2")
spatial = CreateDimReducObject(embeddings = spatial.dim, key = "spatial_")
merfish@reductions[["spatial"]] = spatial

DimPlot(merfish, reduction = "spatial", raster = FALSE)

#QC plots
cols = viridisLite::turbo(n = 256)
p8 = FeaturePlot(merfish, c("nCount_RNA", "nFeature_RNA","volume", "blank_sum", "blank_mean"),
                 reduction = "spatial", pt.size = 0.1, min.cutoff = "q1", max.cutoff = "q99")& scale_colour_gradientn(colors = cols)
p8

FeatureScatter(merfish, "nCount_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "nCount_RNA")
FeatureScatter(merfish, "nFeature_RNA", "blank_mean")

ggplot(merfish@meta.data, aes(log(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(log2(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log2(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(nCount_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nCount_RNA, c(0.01,0.99)))
quantile(merfish$nCount_RNA, c(0.01,0.99))
median(merfish$nCount_RNA)

ggplot(merfish@meta.data, aes(nFeature_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nFeature_RNA, c(0.01,0.99)))
quantile(merfish$nFeature_RNA, c(0.01,0.99))
median(merfish$nFeature_RNA)

ggplot(merfish@meta.data, aes(blank_mean))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$blank_mean, c(0.01,0.99)))
quantile(merfish$blank_mean, c(0.01,0.99))
median(merfish$blank_mean)

#create QC column labeling good and bad cells
good.cells = merfish@meta.data %>% filter(volume <2500 & volume >=40 & nCount_RNA <= 2500 & nCount_RNA >= 30 & nFeature_RNA >= 5 & blank_mean <= 1)
merfish$QCpass = ifelse(colnames(merfish) %in% rownames(good.cells), "TRUE", "FALSE")
table(merfish$QCpass)


## Data subsetting, normalization, scaling------------
merfish = subset(merfish, cells = rownames(good.cells))

merfish <-  NormalizeData(merfish, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish <- ScaleData(merfish, features = rownames(merfish), assay = "RNA", vars.to.regress = "nCount_RNA")

saveRDS(merfish, "data/merfish_b2s20.rds")

#======================Repeat for section b4s8================================

#Read in data-------------------------------------
setwd("~/GOKCE LAB/Code/MERFISH_LPC/github/")
b4s8_matrix = read.csv("data/b4s8/cell_by_gene.csv")
b4s8_matrix[1:10, 1:10]
colnames(b4s8_matrix)
b4s8_matrix$cellID = paste0("b4s8_", b4s8_matrix$X)
length(unique(b4s8_matrix$cellID)) == nrow(b4s8_matrix) #check if all cells have unique name
rownames(b4s8_matrix) = b4s8_matrix$cellID
b4s8_matrix$X = NULL
b4s8_matrix$cellID = NULL

b4s8_meta = read.csv("data/b4s8/cell_metadata.csv")
b4s8_meta[1:20, 1:ncol(b4s8_meta)]
b4s8_meta$cellID = paste0("b4s8_", b4s8_meta$X)
length(unique(b4s8_meta$cellID)) == nrow(b4s8_meta) #check if all cells have unique name
rownames(b4s8_meta) = b4s8_meta$cellID
b4s8_meta$X = NULL

table(rownames(b4s8_matrix) == rownames(b4s8_meta))
table(rownames(b4s8_matrix) %in% rownames(b4s8_meta))
b4s8_matrix = b4s8_matrix[order(match(rownames(b4s8_matrix), rownames(b4s8_meta))),] #reorder rows of matrix to match metadata

b4s8_transcripts = fread("data/b4s8/detected_transcripts.csv", data.table = FALSE) #large file, use data.table::fread() to read it fast
b4s8_transcripts[1:100, 1:10]
b4s8_transcripts$uniqueID = paste0("b4s8_", rep(1:nrow(b4s8_transcripts)))
rownames(b4s8_transcripts) = b4s8_transcripts$uniqueID
b4s8_transcripts$V1 = NULL

#Subsample transcripts for visualization purposes
theme_set(theme_bw())
transcripts_sampled = sample(rownames(b4s8_transcripts), 100000)
transcripts_sampled = b4s8_transcripts[transcripts_sampled,]
p1 = ggplot()+geom_point(data=transcripts_sampled, aes(x = global_x, y = global_y, color = gene), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p1
#ggsave(plot = p1, filename = "outs/b4s8/transcripts_downsampled.pdf", height = 7, width = 9)

#Plotting transcript and cell positions
p2 = ggplot()+geom_point(data=b4s8_meta, aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)+
  theme(legend.position="none")
p2
#ggsave(plot = p2, filename = "outs/b4s8/cells.pdf", height = 7, width = 9)

# Transform (rotate) coordinates of cells, transcripts and images-------------------------
## Rotate cells
embeddings = as.matrix(b4s8_meta[,c("center_x", "center_y")])
ggplot(as.data.frame(embeddings)) + geom_point(aes(x = center_x, y = center_y), size = 0.1, alpha = 0.4)
embeddings.rotated <- as.data.frame(spdep::Rotation(embeddings, 107*pi/180))
colnames(embeddings.rotated) = c("rotated_x", "rotated_y")
ggplot(embeddings.rotated) + geom_point(aes(x = rotated_x, y = rotated_y), size = 0.1, alpha = 0.4)
#ggsave(filename = "outs/b4s8/cells_rotated.pdf", height = 7, width = 9)

head(embeddings.rotated)
head(b4s8_meta)
embeddings.rotated$cellID = rownames(embeddings.rotated)
b4s8_meta = left_join(b4s8_meta, embeddings.rotated, by = "cellID")
rownames(b4s8_meta) = b4s8_meta$cellID

## Rotate ALL transcripts !--------------------------------------------------
transcripts.rotated = as.matrix(b4s8_transcripts[,c("global_x", "global_y")])
transcripts.rotated = as.data.frame(spdep::Rotation(transcripts.rotated, 107*pi/180))
colnames(transcripts.rotated) = c("rotated_x", "rotated_y")
#Add rotated positions to original transcript files
transcripts.rotated$uniqueID = rownames(transcripts.rotated)
b4s8_transcripts = left_join(b4s8_transcripts, transcripts.rotated, by = "uniqueID")
head(b4s8_transcripts)
rownames(b4s8_transcripts) = b4s8_transcripts$uniqueID

#QC and exploratory-------------------------------------
##calculate blank signals
blank_probes = grepl(pattern = "Blank*",colnames(b4s8_matrix))
table(blank_probes)
table(rownames(b4s8_meta) == rownames(b4s8_matrix)) #Must be all TRUE
b4s8_meta$blank_sum = rowSums(b4s8_matrix[blank_probes])
b4s8_meta$blank_mean = rowMeans(b4s8_matrix[blank_probes])

##create Seurat
merfish = CreateSeuratObject(counts = t(b4s8_matrix[!blank_probes]), project = "b4s8", meta.data = b4s8_meta)
p7 = VlnPlot(merfish, features = c("nCount_RNA","nFeature_RNA", "volume", "blank_sum", "blank_mean"), pt.size = 0)
p7
#ggsave(plot = p7, filename = "outs/b4s8/QC_violin.pdf", height = 7, width = 9)

# Add physical spatial dimension to seurat object as dimensional reduction embedding
spatial.dim = as.matrix(merfish@meta.data[c("rotated_x", "rotated_y")])
colnames(spatial.dim) = c("spatial_1", "spatial_2")
spatial = CreateDimReducObject(embeddings = spatial.dim, key = "spatial_")
merfish@reductions[["spatial"]] = spatial

DimPlot(merfish, reduction = "spatial", raster = FALSE)

#QC plots
cols = viridisLite::turbo(n = 256)
p8 = FeaturePlot(merfish, c("nCount_RNA", "nFeature_RNA","volume", "blank_sum", "blank_mean"),
                 reduction = "spatial", pt.size = 0.1, min.cutoff = "q1", max.cutoff = "q99")& scale_colour_gradientn(colors = cols)
p8

FeatureScatter(merfish, "nCount_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "volume")
FeatureScatter(merfish, "nFeature_RNA", "nCount_RNA")
FeatureScatter(merfish, "nFeature_RNA", "blank_mean")

ggplot(merfish@meta.data, aes(log(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(log2(volume)))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(log2(merfish$volume), c(0.01,0.99)))
quantile(merfish$volume, c(0.01,0.99))

ggplot(merfish@meta.data, aes(nCount_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nCount_RNA, c(0.01,0.99)))
quantile(merfish$nCount_RNA, c(0.01,0.99))
median(merfish$nCount_RNA)

ggplot(merfish@meta.data, aes(nFeature_RNA))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$nFeature_RNA, c(0.01,0.99)))
quantile(merfish$nFeature_RNA, c(0.01,0.99))
median(merfish$nFeature_RNA)

ggplot(merfish@meta.data, aes(blank_mean))+stat_ecdf(geom="point")+coord_flip()+
  geom_vline(xintercept = quantile(merfish$blank_mean, c(0.01,0.99)))
quantile(merfish$blank_mean, c(0.01,0.99))
median(merfish$blank_mean)

#create QC column labeling good and bad cells
good.cells = merfish@meta.data %>% filter(volume <2500 & volume >=40 & nCount_RNA <= 2500 & nCount_RNA >= 30 & nFeature_RNA >= 5 & blank_mean <= 1)
merfish$QCpass = ifelse(colnames(merfish) %in% rownames(good.cells), "TRUE", "FALSE")
table(merfish$QCpass)


## Data subsetting, normalization, scaling------------
merfish = subset(merfish, cells = rownames(good.cells))

merfish <-  NormalizeData(merfish, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
merfish <- ScaleData(merfish, features = rownames(merfish), assay = "RNA", vars.to.regress = "nCount_RNA")

saveRDS(merfish, "data/merfish_b4s8.rds")











