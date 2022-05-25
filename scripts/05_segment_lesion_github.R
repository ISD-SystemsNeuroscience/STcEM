library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
setwd(".")

#B3S20----------------------------------------
merfish_b3s20 <- readRDS("data/merfish_b3s20.rds")

df = merfish_b3s20@meta.data %>% select(rotated_x, rotated_y, cell.type_manual)
df = df %>% filter(rotated_x >= -1600 & rotated_x <= 500 & rotated_y >= -5000 & rotated_y <= -3500)

#Read in segmented area polygons
manual.regions2 = readRDS(file = "data/polygons_b3s20.rds")

#Visualize polygons on top of feature or dim plots to verify positioning
# first convert all polygons to one dataframe for simple parsing to ggplot (wont need to add layer for each, can just specify it in one layer)
manual.regions.dfs = list()
for(name in names(manual.regions2)){
  temp =rep(name, nrow(manual.regions2[[name]]))
  manual.regions.dfs[[name]] = cbind(manual.regions2[[name]], temp)
  colnames(manual.regions.dfs[[name]]) = c("rotated_x", "rotated_y", "regionID")
}
manual.regions.dfs = lapply(manual.regions.dfs, as.data.frame)
manual.regions.dfs = bind_rows(manual.regions.dfs) 

#Plot polygons (on any ggplot, or can use Seurat plotting with "combine = F")
#cols3 = c("red","green","blue")
DimPlot(merfish_b3s20, reduction = "spatial", raster = FALSE)

merfish_subset = subset(merfish_b3s20, cells = rownames(df))
p = FeaturePlot(merfish_subset, features = c("Plin2"),
                reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
               fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)
p = DimPlot(merfish_subset, group.by = "cell.cluster_manual",
                reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
                     fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)


#Identify points overlapping polygons, one polygon at a time
library(sp)
poly.list.all = lapply(manual.regions2, Polygon) #convert coordinates to Polygon class
points = SpatialPoints(merfish_b3s20@meta.data %>% select(rotated_x, rotated_y)) #convert cell positions to "SpatialPoints" object
overlap.points = list()
for(region in names(poly.list.all)){
  poly.list = list(poly.list.all[[region]]) 
  poly.list = Polygons(poly.list, "ID") 
  poly.list = list(poly.list) 
  poly.list = SpatialPolygons(poly.list) 
  overlap.points[[region]] = points[poly.list,] #subset points overlapping polygon
}
plot(overlap.points$Core.expanded2) #check if it makes sense

#Label the cells according to the identified overlaps
overlap.cells = list()
for(region in names(overlap.points)){
  overlap.cells[[region]] = rownames(overlap.points[[region]]@coords)
}
#Clean it into non-overlapping regions
lesion.cells = list(lesion.core = overlap.cells$Core,
                          lesion.edge.inner = setdiff(overlap.cells$Core.expanded, overlap.cells$Core),
                          lesion.edge.outer = setdiff(overlap.cells$Core.expanded2, c(overlap.cells$Core, overlap.cells$Core.expanded)),
                          control.WM = overlap.cells$Control.WM,
                          control.GM = setdiff(overlap.cells$Control.GM, overlap.cells$Control.WM))

#Convert regional annotation into dataframes
lesion.cells.dfs = list()
for(name in names(lesion.cells)){
  temp =rep(name, length(lesion.cells[[name]]))
  lesion.cells.dfs[[name]] = as.data.frame(cbind(lesion.cells[[name]], temp))
  colnames(lesion.cells.dfs[[name]]) = c("cellID","region.lesion")
}
lesion.cells.dfs = bind_rows(lesion.cells.dfs)
rownames(lesion.cells.dfs) = lesion.cells.dfs$cellID
lesion.cells.dfs$cellID = NULL
merfish_b3s20 = AddMetaData(merfish_b3s20,metadata = lesion.cells.dfs)
DimPlot(merfish_b3s20, reduction = "spatial", group.by = "region.lesion")
merfish_subset = AddMetaData(merfish_subset,metadata = lesion.cells.dfs)
DimPlot(merfish_subset, reduction = "spatial", group.by = "region.lesion", pt.size = 1)

merfish_b3s20$region.lesion = ifelse(is.na(merfish_b3s20$region.lesion), "Other.area", merfish_b3s20$region.lesion)
saveRDS(merfish_b3s20, file = "data/merfish_b3s20.rds")

#tempDF = left_join(region.fine.cells.dfs, region.coarse.cells.dfs, by = "cellID")

#b2s20----------------------------------------
merfish_b2s20 <- readRDS("data/merfish_b2s20.rds")

df = merfish_b2s20@meta.data %>% select(rotated_x, rotated_y, cell.type_manual)
df = df %>% filter(rotated_x >= -6600 & rotated_x <= -4500 & rotated_y >= -1400 & rotated_y <= 100)

#Read in segmented area polygons
manual.regions2 = readRDS(file = "data/polygons_b2s20.rds")

#Visualize polygons on top of feature or dim plots to verify positioning
# first convert all polygons to one dataframe for simple parsing to ggplot (wont need to add layer for each, can just specify it in one layer)
manual.regions.dfs = list()
for(name in names(manual.regions2)){
  temp =rep(name, nrow(manual.regions2[[name]]))
  manual.regions.dfs[[name]] = cbind(manual.regions2[[name]], temp)
  colnames(manual.regions.dfs[[name]]) = c("rotated_x", "rotated_y", "regionID")
}
manual.regions.dfs = lapply(manual.regions.dfs, as.data.frame)
manual.regions.dfs = bind_rows(manual.regions.dfs) 

#Plot polygons (on any ggplot, or can use Seurat plotting with "combine = F")
#cols3 = c("red","green","blue")
DimPlot(merfish_b2s20, reduction = "spatial", raster = FALSE)

merfish_subset = subset(merfish_b2s20, cells = rownames(df))
p = FeaturePlot(merfish_subset, features = c("Plin2"),
                reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
                     fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)
p = DimPlot(merfish_subset, group.by = "cell.cluster_manual",
            reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
                     fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)


#Identify points overlapping polygons, one polygon at a time
library(sp)
poly.list.all = lapply(manual.regions2, Polygon) #convert coordinates to Polygon class
points = SpatialPoints(merfish_b2s20@meta.data %>% select(rotated_x, rotated_y)) #convert cell positions to "SpatialPoints" object
overlap.points = list()
for(region in names(poly.list.all)){
  poly.list = list(poly.list.all[[region]]) 
  poly.list = Polygons(poly.list, "ID") 
  poly.list = list(poly.list) 
  poly.list = SpatialPolygons(poly.list) 
  overlap.points[[region]] = points[poly.list,] #subset points overlapping polygon
}
plot(overlap.points$Core.expanded2) #check if it makes sense

#Label the cells according to the identified overlaps
overlap.cells = list()
for(region in names(overlap.points)){
  overlap.cells[[region]] = rownames(overlap.points[[region]]@coords)
}
#Clean it into non-overlapping regions
lesion.cells = list(lesion.core = overlap.cells$Core,
                    lesion.edge.inner = setdiff(overlap.cells$Core.expanded, overlap.cells$Core),
                    lesion.edge.outer = setdiff(overlap.cells$Core.expanded2, c(overlap.cells$Core, overlap.cells$Core.expanded)),
                    control.WM = overlap.cells$Control.WM,
                    control.GM = setdiff(overlap.cells$Control.GM, overlap.cells$Control.WM))

#Convert regional annotation into dataframes
lesion.cells.dfs = list()
for(name in names(lesion.cells)){
  temp =rep(name, length(lesion.cells[[name]]))
  lesion.cells.dfs[[name]] = as.data.frame(cbind(lesion.cells[[name]], temp))
  colnames(lesion.cells.dfs[[name]]) = c("cellID","region.lesion")
}
lesion.cells.dfs = bind_rows(lesion.cells.dfs)
rownames(lesion.cells.dfs) = lesion.cells.dfs$cellID
lesion.cells.dfs$cellID = NULL
merfish_b2s20 = AddMetaData(merfish_b2s20,metadata = lesion.cells.dfs)
DimPlot(merfish_b2s20, reduction = "spatial", group.by = "region.lesion")
merfish_subset = AddMetaData(merfish_subset,metadata = lesion.cells.dfs)
DimPlot(merfish_subset, reduction = "spatial", group.by = "region.lesion", pt.size = 1)

merfish_b2s20$region.lesion = ifelse(is.na(merfish_b2s20$region.lesion), "Other.area", merfish_b2s20$region.lesion)
saveRDS(merfish_b2s20, file = "data/merfish_b2s20.rds")


#b4s8----------------------------------------
merfish_b4s8 <- readRDS("data/merfish_b4s8.rds")

df = merfish_b4s8@meta.data %>% select(rotated_x, rotated_y, cell.type_manual)
df = df %>% filter(rotated_x >= -8300 & rotated_x <= -6200 & rotated_y >= 3700 & rotated_y <= 5200)

#Read in segmented area polygons
manual.regions2 = readRDS(file = "data/polygons_b4s8.rds")

#Visualize polygons on top of feature or dim plots to verify positioning
# first convert all polygons to one dataframe for simple parsing to ggplot (wont need to add layer for each, can just specify it in one layer)
manual.regions.dfs = list()
for(name in names(manual.regions2)){
  temp =rep(name, nrow(manual.regions2[[name]]))
  manual.regions.dfs[[name]] = cbind(manual.regions2[[name]], temp)
  colnames(manual.regions.dfs[[name]]) = c("rotated_x", "rotated_y", "regionID")
}
manual.regions.dfs = lapply(manual.regions.dfs, as.data.frame)
manual.regions.dfs = bind_rows(manual.regions.dfs) 

#Plot polygons (on any ggplot, or can use Seurat plotting with "combine = F")
#cols3 = c("red","green","blue")
DimPlot(merfish_b4s8, reduction = "spatial", raster = FALSE)

merfish_subset = subset(merfish_b4s8, cells = rownames(df))
p = FeaturePlot(merfish_subset, features = c("Plin2"),
                reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
                     fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)
p = DimPlot(merfish_subset, group.by = "cell.cluster_manual",
            reduction = "spatial", combine = F, pt.size = 1)
p[[1]]+ geom_polygon(data=manual.regions.dfs, mapping=aes(x=rotated_x, y=rotated_y, group = regionID),
                     fill = NA, alpha = 1, color = "black",linetype="dotted") #+ scale_colour_gradientn(colors = cols3)


#Identify points overlapping polygons, one polygon at a time
library(sp)
poly.list.all = lapply(manual.regions2, Polygon) #convert coordinates to Polygon class
points = SpatialPoints(merfish_b4s8@meta.data %>% select(rotated_x, rotated_y)) #convert cell positions to "SpatialPoints" object
overlap.points = list()
for(region in names(poly.list.all)){
  poly.list = list(poly.list.all[[region]]) 
  poly.list = Polygons(poly.list, "ID") 
  poly.list = list(poly.list) 
  poly.list = SpatialPolygons(poly.list) 
  overlap.points[[region]] = points[poly.list,] #subset points overlapping polygon
}
plot(overlap.points$Core.expanded2) #check if it makes sense

#Label the cells according to the identified overlaps
overlap.cells = list()
for(region in names(overlap.points)){
  overlap.cells[[region]] = rownames(overlap.points[[region]]@coords)
}
#Clean it into non-overlapping regions
lesion.cells = list(lesion.core = overlap.cells$Core,
                    lesion.edge.inner = setdiff(overlap.cells$Core.expanded, overlap.cells$Core),
                    lesion.edge.outer = setdiff(overlap.cells$Core.expanded2, c(overlap.cells$Core, overlap.cells$Core.expanded)),
                    control.WM = overlap.cells$Control.WM,
                    control.GM = setdiff(overlap.cells$Control.GM, overlap.cells$Control.WM))

#Convert regional annotation into dataframes
lesion.cells.dfs = list()
for(name in names(lesion.cells)){
  temp =rep(name, length(lesion.cells[[name]]))
  lesion.cells.dfs[[name]] = as.data.frame(cbind(lesion.cells[[name]], temp))
  colnames(lesion.cells.dfs[[name]]) = c("cellID","region.lesion")
}
lesion.cells.dfs = bind_rows(lesion.cells.dfs)
rownames(lesion.cells.dfs) = lesion.cells.dfs$cellID
lesion.cells.dfs$cellID = NULL
merfish_b4s8 = AddMetaData(merfish_b4s8,metadata = lesion.cells.dfs)
DimPlot(merfish_b4s8, reduction = "spatial", group.by = "region.lesion")
merfish_subset = AddMetaData(merfish_subset,metadata = lesion.cells.dfs)
DimPlot(merfish_subset, reduction = "spatial", group.by = "region.lesion", pt.size = 1)

merfish_b4s8$region.lesion = ifelse(is.na(merfish_b4s8$region.lesion), "Other.area", merfish_b4s8$region.lesion)
saveRDS(merfish_b4s8, file = "data/merfish_b4s8.rds")


#MERGE-----------------
merfish <- readRDS("data/merfish_integrated.rds")

meta_b3s20 = merfish_b3s20@meta.data %>% select("region.lesion")
meta_b2s20 = merfish_b2s20@meta.data %>% select("region.lesion")
meta_b4s8 = merfish_b4s8@meta.data %>% select("region.lesion")
meta = rbind(meta_b3s20, meta_b2s20, meta_b4s8)
merfish = AddMetaData(merfish, metadata = meta)
DimPlot(merfish, group.by = "region.lesion")

saveRDS(merfish, file = "data/merfish_integrated.rds")




